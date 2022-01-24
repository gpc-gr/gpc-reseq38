#
#
#

version 1.0

import "./region2chunks.wdl"


workflow gvcf2vcf_multi_sentieon {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        Array[File] sample_gvcfs
        Array[File] sample_gvcf_indexes

        String region_contig
        Int region_start
        Int region_end
        Int chunk_size = 3 * 1000 * 1000    # 3Mb
        Int chunk_padding = 3000            # 3kb

        String output_prefix

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String python_docker_image = "python:3.8.6-slim-buster"
        String sentieon_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-sentieon:bcftools_1.11-sentieon_202010.02"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1

        String genotype_gvcfs_docker_image = sentieon_docker_image
        Int genotype_gvcfs_threads = 16
        Float genotype_gvcfs_memory_gb = 90

        String vcfconcat_docker_image = sentieon_docker_image
        Int vcfconcat_threads = 4
        Float vcfconcat_memory_gb = 8
    }

    # --------------------------------------------------------------------------------
    # joint genotyping
    # --------------------------------------------------------------------------------

    call region2chunks.region2chunks as step0001_split_region { input:
        contig = region_contig,
        start = region_start,
        end = region_end,
        chunk_size = chunk_size,
        chunks_tsv__name = "${output_prefix}.chunks.tsv",
        docker_image = split_region_docker_image,
        memory_gb = split_region_memory_gb
    }

    scatter (chunk in step0001_split_region.chunks) {

        call step0002_genotype_gvcfs_sentieon { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            individual_gvcfs = sample_gvcfs,
            individual_gvcf_indexes = sample_gvcf_indexes,
            region = "${chunk.contig}:${chunk.start}-${chunk.end}",
            padding = chunk_padding,
            merged_vcf_gz__name = "${output_prefix}.${chunk.name}.vcf.gz",
            docker_image = genotype_gvcfs_docker_image,
            threads = genotype_gvcfs_threads,
            memory_gb = genotype_gvcfs_memory_gb
        }

    }

    call step0003_concat_chunks { input:
        source_vcfs = step0002_genotype_gvcfs_sentieon.merged_vcf_gz,
        source_vcf_indexes = step0002_genotype_gvcfs_sentieon.merged_vcf_gz_index,
        result_vcf_gz__name = "${output_prefix}.vcf.gz",
        docker_image = vcfconcat_docker_image,
        threads = vcfconcat_threads,
        memory_gb = vcfconcat_memory_gb
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File vcf_gz = step0003_concat_chunks.result_vcf_gz
        File vcf_gz_index = step0003_concat_chunks.result_vcf_gz_index
    }

}


task step0002_genotype_gvcfs_sentieon {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        Array[File] individual_gvcfs
        Array[File] individual_gvcf_indexes

        String region
        Int padding

        String merged_vcf_gz__name

        String docker_image
        Int threads
        Float memory_gb
    }

    parameter_meta {
        individual_gvcfs: {
            localization_optional: true
        }
        individual_gvcf_indexes: {
            localization_optional: true
        }
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail
        ulimit -n $(awk '{ print int($0 / 100) }' /proc/sys/fs/file-max)
        ulimit -s 100000

        # step1: extracts target region from individual gVCFs
        python2.7 - <<EOS
        if True:
            import gzip
            import os

            def _get_config_length(template_vcf_gz, target_contig):
                with gzip.open(template_vcf_gz, 'rt') as fin:
                    for line in fin:
                        if not line.startswith('#'):
                            raise Exception

                        if line.startswith('##contig='):
                            kvs = dict(e.split('=', 1) for e in line.strip()[10:-1].split(','))
                            if kvs['ID'] ==  target_contig:
                                return int(kvs['length'])

            #
            with open('~{write_lines(individual_gvcfs)}') as fin:
                original_gvcfs = [l.strip() for l in fin if l.strip()]

            contig, start, end = '~{region}'.replace(':', '-').split('-')
            contig_length = _get_config_length(original_gvcfs[0], contig)

            padding = ~{padding}
            start = max(1, int(start) - padding)
            end = min(int(end) + padding, contig_length)
            region = '{}:{}-{}'.format(contig, start, end)

            #
            with open('./gvcfs.command.txt', 'w') as fout_commands, open('./gvcfs.path.txt', 'w') as fout_gvcfs:
                for original_gvcf in original_gvcfs:
                    extracted_gvcf = os.path.join('./gvcfs.d', os.path.basename(original_gvcf + '.extracted.g.vcf.gz'))
                    commands = [
                        'bcftools view --no-version --regions {} --output-type z --output-file {} {}'.format(region, extracted_gvcf, original_gvcf),
                        'bcftools index --tbi {}'.format(extracted_gvcf)
                    ]

                    print >>fout_commands, ' && '.join(commands)
                    print >>fout_gvcfs, extracted_gvcf
        EOS

        mkdir -p ./gvcfs.d
        cat ./gvcfs.command.txt | xargs -P~{threads} -n1 -I%p sh -c "%p"

        # step2: runs sentieon
        # The following options have been added for GATK 4.1 compatibility:
        #     * --genotype_model multinomial
        #     * --emit_conf 30
        #     * --call_conf 30
        #     * --max_alt_alleles 6
        export VCFCACHE_BLOCKSIZE=4096
        sentieon driver \
            --thread_count ~{threads} \
            --traverse_param 10000/200 \
            --reference ~{reference_fasta} \
            --shard ~{region} \
            --algo GVCFtyper \
            --genotype_model multinomial \
            --emit_conf 30 \
            --call_conf 30 \
            --max_alt_alleles 6 \
            ~{merged_vcf_gz__name} \
            - \
        < ./gvcfs.path.txt

        # step3: removes intermediate files
        rm -rf ./gvcfs.d ./gvcfs.command.txt ./gvcfs.path.txt
    >>>

    output {
        File merged_vcf_gz = "${merged_vcf_gz__name}"
        File merged_vcf_gz_index = "${merged_vcf_gz__name}.tbi"
    }

}


task step0003_concat_chunks {

    input {
        Array[File] source_vcfs
        Array[File] source_vcf_indexes
        String result_vcf_gz__name

        String docker_image
        Int threads
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail
        export VCFCACHE_BLOCKSIZE=4096

        # The following options have been added for GATK 4.1 compatibility:
        #     * --genotype_model multinomial
        #     * --emit_conf 30
        #     * --call_conf 30
        #     * --max_alt_alleles 6
        sentieon driver \
            --thread_count ~{threads} \
            --traverse_param 10000/200 \
            --passthru \
            --algo GVCFtyper \
            --merge \
            --genotype_model multinomial \
            --call_conf 30 \
            --emit_conf 30 \
            --max_alt_alleles 6 \
            ~{result_vcf_gz__name} \
            ~{sep=" " source_vcfs}
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}
