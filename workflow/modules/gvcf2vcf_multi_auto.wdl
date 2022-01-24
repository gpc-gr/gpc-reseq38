#
#
#

version 1.0

import "./gvcf2vcf_multi_cpu.wdl"
import "./gvcf2vcf_multi_sentieon.wdl"


workflow gvcf2vcf_multi_auto {

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

        String? accelaration

        Int genomicsdb_import_batch_size = 256

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"
        String sentieon_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-sentieon:bcftools_1.11-sentieon_202010.02"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1

        String genotype_gvcfs_cpu_docker_image = gatk_docker_image
        Int genotype_gvcfs_cpu_threads = 3
        Float genotype_gvcfs_cpu_memory_gb = 6

        String vcfconcat_cpu_docker_image = bcftools_docker_image
        Int vcfconcat_cpu_threads = 4
        Float vcfconcat_cpu_memory_gb = 1

        String genotype_gvcfs_sentieon_docker_image = sentieon_docker_image
        Int genotype_gvcfs_sentieon_threads = 16
        Float genotype_gvcfs_sentieon_memory_gb = 90

        String vcfconcat_sentieon_docker_image = sentieon_docker_image
        Int vcfconcat_sentieon_threads = 4
        Float vcfconcat_sentieon_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # joint genotyping
    # --------------------------------------------------------------------------------

    if (!defined(accelaration)) {

        call gvcf2vcf_multi_cpu.gvcf2vcf_multi_cpu as step0001_gvcf2vcf_cpu { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            sample_ids = sample_ids,
            sample_gvcfs = sample_gvcfs,
            sample_gvcf_indexes = sample_gvcf_indexes,
            region_contig = region_contig,
            region_start = region_start,
            region_end = region_end,
            chunk_size = chunk_size,
            chunk_padding = chunk_padding,
            output_prefix = output_prefix,
            genomicsdb_import_batch_size = genomicsdb_import_batch_size,
            python_docker_image = python_docker_image,
            split_region_docker_image = split_region_docker_image,
            split_region_memory_gb = split_region_memory_gb,
            genotype_gvcfs_docker_image = genotype_gvcfs_cpu_docker_image,
            genotype_gvcfs_threads = genotype_gvcfs_cpu_threads,
            genotype_gvcfs_memory_gb = genotype_gvcfs_cpu_memory_gb,
            vcfconcat_docker_image = vcfconcat_cpu_docker_image,
            vcfconcat_threads = vcfconcat_cpu_threads,
            vcfconcat_memory_gb = vcfconcat_cpu_memory_gb,
        }

    }

    if (defined(accelaration) && accelaration == "SENTIEON") {

        call gvcf2vcf_multi_sentieon.gvcf2vcf_multi_sentieon as step0001_gvcf2vcf_sentieon { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            sample_gvcfs = sample_gvcfs,
            sample_gvcf_indexes = sample_gvcf_indexes,
            region_contig = region_contig,
            region_start = region_start,
            region_end = region_end,
            chunk_size = chunk_size,
            chunk_padding = chunk_padding,
            output_prefix = output_prefix,
            sentieon_docker_image = sentieon_docker_image,
            split_region_docker_image = split_region_docker_image,
            split_region_memory_gb = split_region_memory_gb,
            genotype_gvcfs_docker_image = genotype_gvcfs_sentieon_docker_image,
            genotype_gvcfs_threads = genotype_gvcfs_sentieon_threads,
            genotype_gvcfs_memory_gb = genotype_gvcfs_sentieon_memory_gb,
            vcfconcat_docker_image = vcfconcat_sentieon_docker_image,
            vcfconcat_threads = vcfconcat_sentieon_threads,
            vcfconcat_memory_gb = vcfconcat_sentieon_memory_gb,
        }

    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File vcf_gz = select_first([
            step0001_gvcf2vcf_cpu.vcf_gz,
            step0001_gvcf2vcf_sentieon.vcf_gz
        ])
        File vcf_gz_index = select_first([
            step0001_gvcf2vcf_cpu.vcf_gz_index,
            step0001_gvcf2vcf_sentieon.vcf_gz_index
        ])
    }

}
