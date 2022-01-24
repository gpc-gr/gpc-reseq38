#
#
#

version 1.0

import "./modules/md5sum.wdl"
import "./modules/region2chunks.wdl"
import "./modules/vcfconcat.wdl"


workflow GPCReseq38_0013_Joint_21_BasicAnnotation {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        String analysis_id

        File reference_fasta
        ResourceVCF dbsnp_vcf
        ResourceVCF clinvar_vcf
        SnpeffDatabase snpeff_database

        Array[String] snpeff_info_keys = ["ANN", "LOF", "NMD"]
        Array[String] clinvar_info_keys = [
            "ALLELEID", "CLNDN", "CLNDNINCL", "CLNDISDBINCL", "CLNHGVS", "CLNREVSTAT", "CLNSIG",
            "CLNSIGCONF", "CLNSIGINCL", "CLNVC", "CLNVCSO", "CLNVI"
        ]

        Array[InputVCF] source_vcfs
        Int chunk_size = 3 * 1000 * 1000   # 3Mb

        Boolean add_af = true
        Boolean add_gf = true
        Array[String] excluded_freq_keys = []

        File sample_idtable = "/dev/null"
        Array[String] sample_category_keys = []
        Array[String] platform_bias_calculation_targets = []
        Float platform_bias_pvalue_threshold = 0.001

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String python_docker_image = "python:3.8.6-slim-buster"
        String r_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-r:dplyr_1.0.7-ggplot2_3.3.5"
        String sequencetoolkit_docker_image = "informationsea/sequencetoolkit:0.2.1"
        String snpeff_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-snpeff:bcftools_1.11-snpeff_4.3.1t"
        String snpsift_docker_image = snpeff_docker_image

        String snpeff_build_docker_image = snpeff_docker_image
        Int snpeff_build_threads = 1
        Float snpeff_build_memory_gb = 8
        String collect_af_keys_docker_image = python_docker_image
        Int collect_af_keys_threads = 1
        Float collect_af_keys_memory_gb = 1
        String hwe_docker_image = bcftools_docker_image
        Int hwe_threads = 4
        Float hwe_memory_gb = 4
        String snpeff_annotate_docker_image = snpeff_docker_image
        Int snpeff_annotate_threads = 4
        Float snpeff_annotate_memory_gb = 6
        String snpsift_dbsnp_docker_image = snpsift_docker_image
        Int snpsift_dbsnp_threads = 4
        Float snpsift_dbsnp_memory_gb = 6
        String add_af_docker_image = sequencetoolkit_docker_image
        Int add_af_threads = 4
        Float add_af_memory_gb = 12
        String normalize_af_vcf_docker_image = bcftools_docker_image
        Int normalize_af_vcf_threads = 4
        Float normalize_af_vcf_memory_gb = 4
        String extract_platform_genotype_counts_docker_image = python_docker_image
        Int extract_platform_genotype_counts_threads = 1
        Float extract_platform_genotype_counts_memory_gb = 4
        String calculate_platform_bias_docker_image = r_docker_image
        Int calculate_platform_bias_threads = 1
        Float calculate_platform_bias_memory_gb = 12
        String create_platform_bias_vcf_docker_image = python_docker_image
        Int create_platform_bias_vcf_threads = 1
        Float create_platform_bias_vcf_memory_gb = 8
        String recompress_platform_bias_vcf_docker_image = bcftools_docker_image
        Int recompress_platform_bias_vcf_threads = 4
        Float recompress_platform_bias_vcf_memory_gb = 4
        String merge_annotations_docker_image = bcftools_docker_image
        Int merge_annotations_threads = 4
        Float merge_annotations_memory_gb = 16
        String collect_annotation_keys_docker_image = python_docker_image
        Int collect_annotation_keys_threads = 1
        Float collect_annotation_keys_memory_gb = 1
        String transfer_annotations_docker_image = bcftools_docker_image
        Int transfer_annotations_threads = 4
        Float transfer_annotations_memory_gb = 16
        String concat_vcfs_docker_image = bcftools_docker_image
        Int concat_vcfs_threads = 8
        Float concat_vcfs_memory_gb = 8
    }

    # --------------------------------------------------------------------------------
    # preparation
    # --------------------------------------------------------------------------------

    call step0001_build_snpeff_database { input:
        reference_fasta = reference_fasta,
        database_id = snpeff_database.id,
        database_gtf = snpeff_database.gtf,
        output_prefix = "${analysis_id}.${snpeff_database.id}.snpeff",
        docker_image = snpeff_build_docker_image,
        threads = snpeff_build_threads,
        memory_gb = snpeff_build_memory_gb
    }

    call step0002_collect_af_keys { input:
        sample_idtable = sample_idtable,
        sample_category_keys = sample_category_keys,
        result_txt__name = "${analysis_id}.af_keys.txt",
        add_af = add_af,
        add_gf = add_gf,
        excluded_freq_keys = excluded_freq_keys,
        docker_image = collect_af_keys_docker_image,
        threads = collect_af_keys_threads,
        memory_gb = collect_af_keys_memory_gb
    }

    # --------------------------------------------------------------------------------
    # annotation
    # --------------------------------------------------------------------------------

    scatter (entry in source_vcfs) {

        call region2chunks.region2chunks as step1000_split_region { input:
            contig = entry.contig,
            start = entry.start,
            end = entry.end,
            chunk_size = chunk_size,
            chunks_tsv__name = "${basename(entry.vcf)}.chunks.tsv"
        }

        scatter (chunk in step1000_split_region.chunks) {

            call step1001_hwe { input:
                source_vcf = entry.vcf,
                source_vcf_index = entry.vcf_index,
                region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                result_vcf_gz__name = sub(basename(entry.vcf), ".vcf.gz$", ".${chunk.name}.hwe.vcf.gz"),
                docker_image = hwe_docker_image,
                threads = hwe_threads,
                memory_gb = hwe_memory_gb
            }

            call step1002_snpeff_annotate { input:
                source_vcf = entry.vcf,
                source_vcf_index = entry.vcf_index,
                region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                snpeff_database_id = snpeff_database.id,
                snpeff_database_config = step0001_build_snpeff_database.config,
                snpeff_database_zip = step0001_build_snpeff_database.zip,
                result_vcf_gz__name = sub(basename(entry.vcf), ".vcf.gz$", ".${chunk.name}.snpeff.vcf.gz"),
                docker_image = snpeff_annotate_docker_image,
                threads = snpeff_annotate_threads,
                memory_gb = snpeff_annotate_memory_gb
            }

            call step1000_snpsift as step1003_snpsift_dbsnp { input:
                source_vcf = entry.vcf,
                source_vcf_index = entry.vcf_index,
                region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                database_vcf = dbsnp_vcf.vcf,
                database_vcf_index = dbsnp_vcf.vcf_index,
                result_vcf_gz__name = sub(basename(entry.vcf), ".vcf.gz$", ".${chunk.name}.dbsnp.vcf.gz"),
                id_only = true,
                info_keys = [],
                docker_image = snpsift_dbsnp_docker_image,
                threads = snpsift_dbsnp_threads,
                memory_gb = snpsift_dbsnp_memory_gb
            }

            call step1000_snpsift as step1004_snpsift_clinvar { input:
                source_vcf = entry.vcf,
                source_vcf_index = entry.vcf_index,
                region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                database_vcf = clinvar_vcf.vcf,
                database_vcf_index = clinvar_vcf.vcf_index,
                result_vcf_gz__name = sub(basename(entry.vcf), ".vcf.gz$", ".${chunk.name}.clinvar.vcf.gz"),
                id_only = false,
                info_keys = clinvar_info_keys,
                docker_image = snpsift_dbsnp_docker_image,
                threads = snpsift_dbsnp_threads,
                memory_gb = snpsift_dbsnp_memory_gb
            }

            if (add_af || add_gf) {

                call step1101_add_af { input:
                    source_vcf = entry.vcf,
                    source_vcf_index = entry.vcf_index,
                    region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                    sample_idtable = sample_idtable,
                    sample_category_keys = sample_category_keys,
                    result_vcf_gz__name = sub(basename(entry.vcf), ".vcf.gz$", ".${chunk.name}.af.vcf.gz"),
                    docker_image = add_af_docker_image,
                    threads = add_af_threads,
                    memory_gb = add_af_memory_gb
                }

                call step1111_normalize_af_vcf { input:
                    source_vcf = step1101_add_af.result_vcf_gz,
                    source_vcf_index = step1101_add_af.result_vcf_gz_index,
                    result_vcf_gz__name = sub(basename(step1101_add_af.result_vcf_gz), ".vcf.gz$", ".normalized.vcf.gz"),
                    docker_image = normalize_af_vcf_docker_image,
                    threads = normalize_af_vcf_threads,
                    memory_gb = normalize_af_vcf_memory_gb
                }

                call step1112_extract_platform_genotype_counts { input:
                    vcf = step1111_normalize_af_vcf.result_vcf_gz,
                    vcf_index = step1111_normalize_af_vcf.result_vcf_gz_index,
                    platforms = platform_bias_calculation_targets,
                    tsv__name = "${basename(step1111_normalize_af_vcf.result_vcf_gz)}.platform_genotype_counts.tsv",
                    docker_image = extract_platform_genotype_counts_docker_image,
                    threads = extract_platform_genotype_counts_threads,
                    memory_gb = extract_platform_genotype_counts_memory_gb
                }

                call step1113_calculate_platform_bias { input:
                    source_tsv = step1112_extract_platform_genotype_counts.tsv,
                    platforms = platform_bias_calculation_targets,
                    result_tsv__name = "${basename(step1112_extract_platform_genotype_counts.tsv)}.pvalue.tsv",
                    docker_image = calculate_platform_bias_docker_image,
                    threads = calculate_platform_bias_threads,
                    memory_gb = calculate_platform_bias_memory_gb
                }

                call step1114_create_platform_bias_vcf { input:
                    af_vcf = step1101_add_af.result_vcf_gz,
                    platform_bias_pvalue_tsv = step1113_calculate_platform_bias.result_tsv,
                    platform_bias_pvalue_threshold = platform_bias_pvalue_threshold,
                    result_vcf_gz__name = "${basename(step1101_add_af.result_vcf_gz)}.platform_bias.vcf.gz",
                    docker_image = create_platform_bias_vcf_docker_image,
                    threads = create_platform_bias_vcf_threads,
                    memory_gb = create_platform_bias_vcf_memory_gb
                }

                call step1115_recompress_platform_bias_vcf { input:
                    source_vcf = step1114_create_platform_bias_vcf.result_vcf_gz,
                    result_vcf_gz__name = "${basename(step1101_add_af.result_vcf_gz)}.platform_bias.vcf.gz",
                    docker_image = recompress_platform_bias_vcf_docker_image,
                    threads = recompress_platform_bias_vcf_threads,
                    memory_gb = recompress_platform_bias_vcf_memory_gb
                }

            }

            call step1201_merge_annotations { input:
                dbsnp_vcf = step1003_snpsift_dbsnp.result_vcf_gz,
                dbsnp_vcf_index = step1003_snpsift_dbsnp.result_vcf_gz_index,
                other_source_vcfs = select_all([
                    step1101_add_af.result_vcf_gz,
                    step1115_recompress_platform_bias_vcf.result_vcf_gz,
                    step1001_hwe.result_vcf_gz,
                    step1002_snpeff_annotate.result_vcf_gz,
                    step1004_snpsift_clinvar.result_vcf_gz
                ]),
                other_source_vcf_indexes = select_all([
                    step1101_add_af.result_vcf_gz_index,
                    step1115_recompress_platform_bias_vcf.result_vcf_gz_index,
                    step1001_hwe.result_vcf_gz_index,
                    step1002_snpeff_annotate.result_vcf_gz_index,
                    step1004_snpsift_clinvar.result_vcf_gz_index
                ]),
                result_vcf_gz__name = sub(basename(entry.vcf), ".vcf(.gz)?$", ".${chunk.name}.annotations.vcf.gz"),
                docker_image = merge_annotations_docker_image,
                threads = merge_annotations_threads,
                memory_gb = merge_annotations_memory_gb
            }

            if (step1201_merge_annotations.result_vcf_gz_is_not_empty) {

                Array[String] af_info_keys = if (add_af) then step0002_collect_af_keys.info_keys else []
                Array[String] platform_bias_info_keys = if (add_af) then select_first([step1114_create_platform_bias_vcf.info_keys]) else []

                call step1202_transfer_annotations { input:
                    source_vcf = entry.vcf,
                    source_vcf_index = entry.vcf_index,
                    annotation_vcf = step1201_merge_annotations.result_vcf_gz,
                    annotation_vcf_index = step1201_merge_annotations.result_vcf_gz_index,
                    result_vcf_gz__name = sub(basename(entry.vcf), ".vcf.gz$", ".${chunk.name}.annotated.vcf.gz"),
                    region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                    info_keys = flatten([
                        af_info_keys,
                        platform_bias_info_keys,
                        ["HWE"],
                        snpeff_info_keys,
                        clinvar_info_keys
                    ]),
                    docker_image = transfer_annotations_docker_image,
                    threads = transfer_annotations_threads,
                    memory_gb = transfer_annotations_memory_gb
                }

            }

        }

        call vcfconcat.vcfconcat as step2001_concat_vcfs { input:
            source_vcfs = select_all(step1202_transfer_annotations.result_vcf_gz),
            source_vcf_indexes = select_all(step1202_transfer_annotations.result_vcf_gz_index),
            result_vcf_gz__name = sub(basename(entry.vcf), ".vcf(.gz)?$", ".${dbsnp_vcf.id}_${clinvar_vcf.id}_${snpeff_database.id}.vcf.gz"),
            apply_sort = false,
            source_vcfs_have_same_header = true,
            create_index = true,
            docker_image = concat_vcfs_docker_image,
            threads = concat_vcfs_threads,
            memory_gb = concat_vcfs_memory_gb
        }

    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            step2001_concat_vcfs.result_vcf_gz,
            step2001_concat_vcfs.result_vcf_gz_index
        ]),
        md5sum_txt__name = "${analysis_id}.GPCReseq38_0013_Joint_21_BasicAnnotation.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        Array[File] vcf_gz = step2001_concat_vcfs.result_vcf_gz
        Array[File] vcf_gz_index = step2001_concat_vcfs.result_vcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct ResourceVCF {
    String id
    File vcf
    File vcf_index
}


struct SnpeffDatabase {
    String id
    File gtf
}


struct InputVCF {
    File vcf
    File vcf_index
    String contig
    Int start
    Int end
}


task step0001_build_snpeff_database {

    input {
        File reference_fasta
        String database_id
        File database_gtf
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
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        mkdir -p data/~{database_id}
        ln -s ~{reference_fasta} data/~{database_id}/sequences.fa
        ln -s ~{database_gtf} data/~{database_id}/genes.gtf

        {
            cat /opt/conda/share/snpeff-*/snpEff.config
            echo
            echo '~{database_id}.genome : ~{basename(database_gtf)}'
            echo '~{database_id}.M.codonTable : Vertebrate_Mitochondrial'
            echo '~{database_id}.MT.codonTable : Vertebrate_Mitochondrial'
        } > ~{output_prefix}.config

        snpEff build -c ~{output_prefix}.config -verbose -gtf22 ~{database_id}

        rm data/~{database_id}/sequences.fa
        rm data/~{database_id}/genes.gtf

        zip -0 -r ~{output_prefix}.zip data
        rm -rf data
    >>>

    output {
        File config = "${output_prefix}.config"
        File zip = "${output_prefix}.zip"
    }

}


task step0002_collect_af_keys {

    input {
        File sample_idtable
        Array[String] sample_category_keys
        String result_txt__name

        Boolean add_af
        Boolean add_gf
        Array[String] excluded_freq_keys

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

        python3 - <<EOS >~{result_txt__name}
        if True:
            import collections
            import csv

            #
            with open('~{write_lines(sample_category_keys)}') as fin:
                keys = [k.strip() for k in fin if k.strip()]

            #
            values = collections.defaultdict(set)
            with open('~{sample_idtable}') as fin:
                for record in csv.DictReader(fin, delimiter='\t'):
                    for key in keys:
                        values[key].add(record[key])

            #
            prefixes = []
            if ~{true="True" false="False" add_af}:
                prefixes.extend(['AC', 'AN', 'AF'])
            if ~{true="True" false="False" add_gf}:
                prefixes.extend(['GenotypeCount', 'GenotypeCount', 'nhomalt'])

            info_keys = list(prefixes)
            for key in keys:
                for value in sorted(values[key]):
                    info_keys.extend(f'{p}_{value}' for p in prefixes)

            #
            with open('~{write_lines(excluded_freq_keys)}') as fin:
                excluded_info_keys = set(l.strip() for l in fin if l.strip())

            for info_key in info_keys:
                if info_key not in excluded_info_keys:
                    print(info_key)
        EOS
    >>>

    output {
        File result_txt = "${result_txt__name}"
        Array[String] info_keys = read_lines(result_txt__name)
    }
}


task step1001_hwe {

    input {
        File source_vcf
        File source_vcf_index
        String region
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

        bcftools +fill-tags \
            --no-version \
            --threads ~{threads} \
            --regions ~{region} \
            --targets ~{region} \
            --output-type u \
            -- \
            --tags HWE \
            ~{source_vcf} \
        | bcftools view \
            --no-version \
            --threads ~{threads} \
            --drop-genotypes \
            --output-type u \
        | bcftools annotate \
            --no-version \
            --remove ID,QUAL,FILTER,^INFO/HWE \
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


task step1002_snpeff_annotate {

    input {
        File source_vcf
        File source_vcf_index
        String region
        String snpeff_database_id
        File snpeff_database_config
        File snpeff_database_zip
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
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 6/10)}m"

        #
        cat ~{snpeff_database_config} > ~{basename(snpeff_database_config)}
        unzip ~{snpeff_database_zip}

        #
        bcftools view \
            --no-version \
            --threads ~{threads} \
            --regions ~{region} \
            --targets ~{region} \
            --drop-genotypes \
            --output-type u \
            ~{source_vcf} \
        | bcftools annotate \
            --no-version \
            --threads ~{threads} \
            --remove ID,QUAL,FILTER,^INFO/ \
            --output-type v \
        | snpEff ann \
            -config ~{basename(snpeff_database_config)} \
            -verbose \
            ~{snpeff_database_id} \
        | bcftools view \
            --no-version \
            --threads ~{threads} \
            --output-type z \
            --output-file ~{result_vcf_gz__name}

        bcftools index \
            --threads ~{threads} \
            --tbi \
            ~{result_vcf_gz__name}

        #
        rm -rf ~{basename(snpeff_database_config)}
        rm -rf data
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}


task step1000_snpsift {

    input {
        File source_vcf
        File source_vcf_index
        String region
        File database_vcf
        File database_vcf_index
        String result_vcf_gz__name

        Boolean id_only
        Array[String] info_keys

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
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 6/10)}m"

        #
        bcftools view \
            --no-version \
            --threads ~{threads} \
            --regions ~{region} \
            --targets ~{region} \
            --drop-genotypes \
            --output-type u \
            ~{source_vcf} \
        | bcftools annotate \
            --no-version \
            --remove ID,QUAL,FILTER,^INFO/ \
            --output-type v \
        | SnpSift annotate \
            -v \
            ~{true="-id" false="" id_only} \
            -info '~{sep="," info_keys}' \
            -tabix \
            ~{database_vcf} \
        | bcftools view \
            --no-version \
            --threads ~{threads} \
            --output-type z \
            --output-file ~{result_vcf_gz__name}

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


task step1101_add_af {

    input {
        File source_vcf
        File source_vcf_index
        String region
        File sample_idtable
        Array[String] sample_category_keys
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

        bcftools view \
            --no-version \
            --threads ~{threads} \
            --regions ~{region} \
            --targets ~{region} \
            --output-type u \
            ~{source_vcf} \
        | bcftools annotate \
            --no-version \
            --remove ID,QUAL,FILTER,^INFO/ \
            --output-type v \
        | sequencetoolkit vcfutils add-af \
            --category ~{sample_idtable} \
            --id id \
            --value ~{sep=" " sample_category_keys} \
        | bcftools view \
            --no-version \
            --threads ~{threads} \
            --drop-genotypes \
            --output-type z \
            --output-file ~{result_vcf_gz__name}

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


task step1111_normalize_af_vcf {

    input {
        File source_vcf
        File source_vcf_index
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

        bcftools norm \
            --no-version \
            --threads ~{threads} \
            --multiallelics -any \
            --do-not-normalize \
            --output-type z \
            --output ~{result_vcf_gz__name} \
            ~{source_vcf}

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


task step1112_extract_platform_genotype_counts {

    input {
        File vcf
        File vcf_index
        Array[String] platforms
        String tsv__name

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

        python3 - <<EOS
        if True:
            import gzip

            with open('~{write_lines(platforms)}') as fin:
                PLATFORMS = [l.strip() for l in fin if l.strip()]

            def _parse_vcf(fin):
                for line in fin:
                    if not line.startswith('#'):
                        yield _parse_vcf_line(line)

            def _parse_vcf_line(line):
                cols = line.split('\t')
                contig = cols[0]
                position = cols[1]
                reference = cols[3]
                alternative = cols[4]
                af = [float(e[3:]) for e in cols[7].split(';') if e.startswith('AF=')][0]
                genotype_counts = [
                    [e.split('=')[1] for e in cols[7].split(';') if e.startswith(f'GenotypeCount_{platform}')][0]
                    for platform in PLATFORMS
                ]

                return contig, position, reference, alternative, af, ','.join(genotype_counts)

            def _read_vcf(fin):
                current_contig_and_position = None, None
                current_records = []

                for record in _parse_vcf(fin):
                    if tuple(record[:2]) != current_contig_and_position:
                        if current_records:
                            yield current_records

                        current_contig_and_position = tuple(record[:2])
                        current_records = [record]

                    else:
                        current_records.append(record)

                if current_records:
                    yield current_records

            def _open(path):
                return gzip.open(path, 'rt') if path.endswith('gz') else open(path)

            #
            print('\t'.join(['contig', 'position', 'reference', 'alternative', 'allele_frequency', 'genotype_counts']))

            with _open('~{vcf}') as fin:
                for records in _read_vcf(fin):
                    max_af_record = list(sorted(records, key=lambda r: (-r[4], r[3])))[0]
                    print('\t'.join(map(str, max_af_record)))
        EOS
    >>>

    output {
        File tsv = "${tsv__name}"
    }

}


task step1113_calculate_platform_bias {

    input {
        File source_tsv
        Array[String] platforms
        String result_tsv__name

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

        Rscript - <<EOS
            d <- read.table('~{source_tsv}', head=T, sep='\t')
            d\$pvalue <- apply(d, 1, function (r) {
                m <- matrix(sapply(strsplit(r[6], ',')[[1]], as.numeric), nrow=~{length(platforms)}, byrow=T)
                ret <- try((function () {
                    fisher.test(m, workspace=100000000, simulate.p.value=T)$p.value
                })(), silent=F)

                if (class(ret) == 'try-error') {
                    ret <- NA
                }

                ret
            })

            write.table(d, '~{result_tsv__name}', sep='\t', quote=F, row.names=F)
        EOS
    >>>

    output {
        File result_tsv = "${result_tsv__name}"
    }

}


task step1114_create_platform_bias_vcf {

    input {
        File af_vcf
        File platform_bias_pvalue_tsv
        Float platform_bias_pvalue_threshold
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

        python3 - <<EOS | gzip -c > ~{result_vcf_gz__name}
        if True:
            import csv
            import gzip

            with open('~{platform_bias_pvalue_tsv}') as fin:
                sites = {
                    (r['contig'], r['position']): r
                    for r in csv.DictReader(fin, delimiter='\t')
                }

            #
            platform_bias_pvalue_threshold = ~{platform_bias_pvalue_threshold}
            info_written = False

            for line in gzip.open('~{af_vcf}', 'rt'):
                line = line.strip()
                if not line:
                    continue

                if line.startswith('#'):
                    if line.startswith('##INFO=') and (not info_written):
                        print(f'##INFO=<ID=TOMMO_PLATFORM_BIAS_TEST_ALLELE,Number=1,Type=String,Description="an alternative allele used for platform-bias test">')
                        print(f'##INFO=<ID=TOMMO_PLATFORM_BIAS_TEST_PVALUE,Number=1,Type=Float,Description="p-value of platform-bias test">')
                        print(f'##INFO=<ID=TOMMO_POSSIBLE_PLATFORM_BIAS,Number=0,Type=Flag,Description="p-value of platform-bias test is equal or less than {platform_bias_pvalue_threshold}">')
                        info_written = True

                    print(line)

                else:
                    #
                    cols = line.strip().split('\t')
                    contig = cols[0]
                    position = cols[1]
                    reference = cols[3]
                    alternatives = cols[4].split(',')

                    if (contig, position) not in sites:
                        raise Exception(f'{contig}:{position}')

                    site = sites[contig, position]
                    if (site['alternative'] not in alternatives):
                        raise Exception(f'{contig}:{position}')

                    #
                    info = ['TOMMO_PLATFORM_BIAS_TEST_ALLELE={}'.format(site['alternative'])]
                    if site['pvalue'] == 'NA':
                        info.append('TOMMO_PLATFORM_BIAS_TEST_PVALUE=NA')
                    else:
                        pvalue = float(site['pvalue'])
                        info.append('TOMMO_PLATFORM_BIAS_TEST_PVALUE={:.06f}'.format(pvalue))
                        if pvalue <= platform_bias_pvalue_threshold:
                            info.append('TOMMO_POSSIBLE_PLATFORM_BIAS')

                    #
                    cols[7] = ';'.join(info)
                    print('\t'.join(cols))
        EOS
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        Array[String] info_keys = [
            "TOMMO_PLATFORM_BIAS_TEST_ALLELE",
            "TOMMO_PLATFORM_BIAS_TEST_PVALUE",
            "TOMMO_POSSIBLE_PLATFORM_BIAS"
        ]
    }

}


task step1115_recompress_platform_bias_vcf {

    input {
        File source_vcf
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

        bcftools view \
            --no-version \
            --threads ~{threads} \
            --output-type z \
            --output-file ~{result_vcf_gz__name} \
            ~{source_vcf}

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


task step1201_merge_annotations {

    input {
        File dbsnp_vcf
        File dbsnp_vcf_index
        Array[File] other_source_vcfs
        Array[File] other_source_vcf_indexes
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

        #
        bcftools merge \
            --no-version \
            --threads ~{threads} \
            --file-list ~{write_lines(other_source_vcfs)} \
            --output-type b \
            --output ~{result_vcf_gz__name}.temp.bcf.gz

        bcftools index \
            --threads ~{threads} \
            --csi \
            ~{result_vcf_gz__name}.temp.bcf.gz

        #
        bcftools annotate \
            --no-version \
            --threads ~{threads} \
            --annotations ~{dbsnp_vcf} \
            --columns ID \
            --output-type z \
            --output ~{result_vcf_gz__name} \
            ~{result_vcf_gz__name}.temp.bcf.gz

        bcftools index \
            --threads ~{threads} \
            --tbi \
            ~{result_vcf_gz__name}

        #
        rm ~{result_vcf_gz__name}.temp.bcf.gz*

        #
        num_records=$(bcftools view --no-header ~{result_vcf_gz__name} | head -n1 | wc -l || true)
        if [ ${num_records} -eq 1 ]; then
            echo true > ~{result_vcf_gz__name}.is_not_empty.txt
        else
            echo false > ~{result_vcf_gz__name}.is_not_empty.txt
        fi
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
        Boolean result_vcf_gz_is_not_empty = read_boolean("${result_vcf_gz__name}.is_not_empty.txt")
    }

}


task step1202_collect_annotation_keys {

    input {
        File annotation_vcf
        File annotation_vcf_index

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

        python3 - <<EOS >info.txt
        if True:
            import gzip
            import re

            with gzip.open('~{annotation_vcf}', 'rt') as fin:
                for line in fin:
                    if not line.startswith('#'):
                        break

                    match = re.search('##INFO=<ID=([a-zA-Z0-9]+),>')
                    if match:
                        print(match.group(1))
        EOS
    >>>

    output {
        Array[String] info_keys = read_lines("info.txt")
    }

}


task step1202_transfer_annotations {

    input {
        File source_vcf
        File source_vcf_index
        File annotation_vcf
        File annotation_vcf_index
        String result_vcf_gz__name

        String region
        Array[String] info_keys

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

        #
        bcftools view \
            --no-version \
            --threads ~{threads} \
            --regions ~{region} \
            --targets ~{region} \
            --output-type b \
            --output-file ~{result_vcf_gz__name}.target.bcf.gz \
            ~{source_vcf}

        bcftools index \
            --threads ~{threads} \
            --csi \
            ~{result_vcf_gz__name}.target.bcf.gz

        #
        bcftools annotate \
            --no-version \
            --threads ~{threads} \
            --annotations ~{annotation_vcf} \
            --columns ID,~{sep="," prefix("INFO/", info_keys)} \
            --output-type z \
            --output ~{result_vcf_gz__name} \
            ~{result_vcf_gz__name}.target.bcf.gz

        bcftools index \
            --threads ~{threads} \
            --tbi \
            ~{result_vcf_gz__name}

        #
        rm ~{result_vcf_gz__name}.target.bcf.gz*
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}
