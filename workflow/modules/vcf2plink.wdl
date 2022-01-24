#
#
#

version 1.0

import "./region2chunks.wdl"
import "./vcfconcat.wdl"


workflow vcf2plink {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        Array[File] source_vcfs
        Array[File] source_vcf_indexes
        Array[String] contigs
        Array[Int] starts
        Array[Int] ends
        String output_prefix

        File? target_sample_list

        Boolean pass_only = false
        Boolean snv_only = false
        Boolean biallelic_only = false

        Float? hwe_ge
        Float? maf_ge
        Float? geno_le

        Boolean normalize_position_and_alleles = true

        Int chunk_size = 3 * 1000 * 1000

        File? sample_idtable
        String? sample_sex_column_key

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String plink_docker_image = "quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"
        String python_docker_image = "python:3.8.6-slim-buster"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1
        String make_slim_vcf_docker_image = bcftools_docker_image
        Int make_slim_vcf_threads = 4
        Float make_slim_vcf_memory_gb = 16
        String concat_vcfs_docker_image = bcftools_docker_image
        Int concat_vcfs_threads = 8
        Float concat_vcfs_memory_gb = 16
        String vcf2plink_docker_image = plink_docker_image
        Int vcf2plink_threads = 1
        Float vcf2plink_memory_gb = 64
        String update_sex_in_fam_docker_image = python_docker_image
        Float update_sex_in_fam_memory_gb = 4
    }

    scatter (index in range(length(source_vcfs))) {

        call region2chunks.region2chunks as step0000_split_region { input:
            contig = contigs[index],
            start = starts[index],
            end = ends[index],
            chunk_size = chunk_size,
            chunks_tsv__name = "${basename(source_vcfs[index])}.chunks.tsv",
            docker_image = split_region_docker_image,
            memory_gb = split_region_memory_gb
        }

        scatter (chunk in step0000_split_region.chunks) {

            call step0001_make_slim_vcf { input:
                reference_fasta = reference_fasta,
                reference_fasta_general_indexes = reference_fasta_general_indexes,
                source_vcf = source_vcfs[index],
                source_vcf_index = source_vcf_indexes[index],
                region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                result_vcf_gz__name = sub(basename(source_vcfs[index]), ".vcf(.gz)?$", ".${chunk.name}.slim.vcf.gz"),
                target_sample_list = target_sample_list,
                pass_only = pass_only,
                snv_only = snv_only,
                biallelic_only = biallelic_only,
                hwe_ge = hwe_ge,
                maf_ge = maf_ge,
                geno_le = geno_le,
                normalize_position_and_alleles = normalize_position_and_alleles,
                docker_image = make_slim_vcf_docker_image,
                threads = make_slim_vcf_threads,
                memory_gb = make_slim_vcf_memory_gb,
            }

        }

    }

    call vcfconcat.vcfconcat as step1001_concat_vcfs { input:
        source_vcfs = flatten(step0001_make_slim_vcf.result_vcf_gz),
        source_vcf_indexes = flatten(step0001_make_slim_vcf.result_vcf_gz_index),
        result_vcf_gz__name = "${output_prefix}.vcf.gz",
        apply_sort = false,
        source_vcfs_have_same_header = true,
        create_index = false,
        docker_image = concat_vcfs_docker_image,
        threads = concat_vcfs_threads,
        memory_gb = concat_vcfs_memory_gb
    }

    call step1002_vcf2plink { input:
        vcf = step1001_concat_vcfs.result_vcf_gz,
        output_prefix = output_prefix,
        docker_image = vcf2plink_docker_image,
        threads = vcf2plink_threads,
        memory_gb = vcf2plink_memory_gb
    }

    if (defined(sample_idtable) && defined(sample_sex_column_key)) {

        call step1101_update_sample_sex_in_fam { input:
            sample_idtable = select_first([sample_idtable]),
            sample_sex_column_key = select_first([sample_sex_column_key]),
            original_fam = step1002_vcf2plink.fam,
            docker_image = update_sex_in_fam_docker_image,
            memory_gb = update_sex_in_fam_memory_gb
        }

    }

    output {
        File bed = step1002_vcf2plink.bed
        File bim = step1002_vcf2plink.bim
        File fam = select_first([step1101_update_sample_sex_in_fam.result_fam, step1002_vcf2plink.fam])
    }

}


task step0001_make_slim_vcf {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        File source_vcf
        File source_vcf_index
        String region
        String result_vcf_gz__name

        File? target_sample_list

        Boolean pass_only
        Boolean snv_only
        Boolean biallelic_only
        Float? hwe_ge
        Float? maf_ge
        Float? geno_le
        Boolean normalize_position_and_alleles

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

        bcftools norm \
            --no-version \
            --threads ~{threads} \
            --regions ~{region} \
            --targets ~{region} \
            --fasta-ref ~{reference_fasta} \
            --multiallelics -any \
            ~{true="" false="--do-not-normalize" normalize_position_and_alleles} \
            --output-type u \
            ~{source_vcf} \
        | bcftools view \
            --no-version \
            ~{"--samples-file " + target_sample_list} \
            ~{true="--apply-filters PASS" false="" pass_only} \
            ~{true="--types snps" false="" snv_only} \
            ~{true="--min-alleles 2" false="" biallelic_only} \
            ~{true="--max-alleles 2" false="" biallelic_only} \
            --output-type u \
        | bcftools +fill-tags \
            --no-version \
            --output-type u \
            -- \
            --tags INFO/HWE,INFO/MAF,INFO/F_MISSING \
        | bcftools filter \
            --no-version \
            --include '1 == 1 ~{"&& INFO/HWE >= " + hwe_ge} ~{"&& INFO/MAF >= " + maf_ge} ~{"&& INFO/F_MISSING <= " + geno_le}' \
            --output-type u \
        | bcftools annotate \
            --no-version \
            --threads ~{threads} \
            --set-id '%CHROM\:%POS\_%REF\_%FIRST_ALT' \
            --remove QUAL,FILTER,^INFO/,^FORMAT/GT \
            --output-type z \
            --output ~{result_vcf_gz__name}

        bcftools index \
            --threads ~{threads} \
            --tbi \
            ~{result_vcf_gz__name}
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}


task step1002_vcf2plink {

    input {
        File vcf
        String output_prefix

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

        plink \
            --vcf ~{vcf} \
            --keep-allele-order \
            --make-bed \
            --out ~{output_prefix} \
            --memory ~{floor(memory_gb * 1000 * 8/10)}
    >>>

    output {
        File bed = "${output_prefix}.bed"
        File bim = "${output_prefix}.bim"
        File fam = "${output_prefix}.fam"
    }

}


task step1101_update_sample_sex_in_fam {

    input {
        File sample_idtable
        String sample_sex_column_key
        File original_fam

        String docker_image
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        python3 - <<EOS > ~{basename(original_fam)}
        if True:
            import csv

            with open('~{sample_idtable}') as fin:
                sample_sex_map = {r['id']: r['~{sample_sex_column_key}'] for r in csv.DictReader(fin, delimiter='\t')}

            with open('~{original_fam}') as fin:
                for line in fin:
                    cols = line.strip().split()
                    if cols[1] in sample_sex_map:
                        cols[4] = sample_sex_map[cols[1]]

                    print(' '.join(cols))
        EOS
    >>>

    output {
        File result_fam = "${basename(original_fam)}"
    }

}
