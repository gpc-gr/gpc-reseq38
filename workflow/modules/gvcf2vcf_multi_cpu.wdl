#
#
#

version 1.0

import "./region2chunks.wdl"
import "./vcfconcat.wdl"


workflow gvcf2vcf_multi_cpu {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        Array[String] sample_ids
        Array[File] sample_gvcfs
        Array[File] sample_gvcf_indexes

        String region_contig
        Int region_start
        Int region_end
        Int chunk_size = 3 * 1000 * 1000    # 3Mb
        Int chunk_padding = 3000            # 3kb

        String output_prefix

        Int genomicsdb_import_batch_size = 256

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1

        String genotype_gvcfs_docker_image = gatk_docker_image
        Int genotype_gvcfs_threads = 3
        Float genotype_gvcfs_memory_gb = 6

        String vcfconcat_docker_image = bcftools_docker_image
        Int vcfconcat_threads = 4
        Float vcfconcat_memory_gb = 1
    }

    # --------------------------------------------------------------------------------
    # joint genotyping
    # --------------------------------------------------------------------------------

    call region2chunks.region2chunks as  step0001_split_region { input:
        contig = region_contig,
        start = region_start,
        end = region_end,
        chunk_size = chunk_size,
        chunks_tsv__name = "${output_prefix}.chunks.tsv",
        docker_image = split_region_docker_image,
        memory_gb = split_region_memory_gb
    }

    scatter (chunk in step0001_split_region.chunks) {

        call step0002_genotype_gvcfs_cpu { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            individual_ids = sample_ids,
            individual_gvcfs = sample_gvcfs,
            individual_gvcf_indexes = sample_gvcf_indexes,
            region = "${chunk.contig}:${chunk.start}-${chunk.end}",
            padding = chunk_padding,
            merged_vcf_gz__name = "${output_prefix}.${chunk.name}.vcf.gz",
            genomicsdb_import_batch_size = genomicsdb_import_batch_size,
            docker_image = genotype_gvcfs_docker_image,
            threads = genotype_gvcfs_threads,
            memory_gb = genotype_gvcfs_memory_gb
        }

    }

    call vcfconcat.vcfconcat as step1001_concat_chunks { input:
        source_vcfs = step0002_genotype_gvcfs_cpu.merged_vcf_gz,
        source_vcf_indexes = step0002_genotype_gvcfs_cpu.merged_vcf_gz_index,
        result_vcf_gz__name = "${output_prefix}.vcf.gz",
        apply_sort = false,
        source_vcfs_have_same_header = true,
        create_index = true,
        docker_image = vcfconcat_docker_image,
        threads = vcfconcat_threads,
        memory_gb = vcfconcat_memory_gb
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File vcf_gz = step1001_concat_chunks.result_vcf_gz
        File vcf_gz_index = step1001_concat_chunks.result_vcf_gz_index
    }

}


task step0002_genotype_gvcfs_cpu {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        Array[String] individual_ids
        Array[File] individual_gvcfs
        Array[File] individual_gvcf_indexes

        String region
        Int padding
        String merged_vcf_gz__name

        Int genomicsdb_import_batch_size

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

        #
        paste ~{write_lines(individual_ids)} ~{write_lines(individual_gvcfs)} > sources.txt

        #
        export JAVA_TOOL_OPTIONS="-XX:+UseSerialGC -XX:CICompilerCount=2 -Xmx~{floor(memory_gb * 1000 * 8/10)}m"
        gatk GenomicsDBImport \
            --reference ~{reference_fasta} \
            --intervals ~{region} \
            --interval-padding ~{padding} \
            --sample-name-map sources.txt \
            --genomicsdb-workspace-path ./genomicsdb-workspace \
            --batch-size ~{genomicsdb_import_batch_size} \
            --reader-threads ~{threads}

        #
        gatk GenotypeGVCFs \
            --reference ~{reference_fasta} \
            --intervals ~{region} \
            --variant gendb://./genomicsdb-workspace \
            --only-output-calls-starting-in-intervals \
            --output ~{merged_vcf_gz__name} \
            --create-output-variant-index

        #
        rm -rf sources.txt
        rm -rf ./genomicsdb-workspace
    >>>

    output {
        File merged_vcf_gz = "${merged_vcf_gz__name}"
        File merged_vcf_gz_index = "${merged_vcf_gz__name}.tbi"
    }

}
