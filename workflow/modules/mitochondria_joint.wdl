#
#
#

version 1.0

import "./gvcf2vcf_multi_auto.wdl"
import "./gvcf2vcf_multi_sentieon.wdl"
import "./mitochondria_merge.wdl"


workflow mitochondria_joint {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File baseline_reference_fasta
        Array[File] baseline_reference_fasta_general_indexes = [
            "${baseline_reference_fasta}.fai",
            sub(baseline_reference_fasta, ".fa(sta)?$", ".dict")
        ]

        File mitochondria_shifted_reference_fasta
        Array[File] mitochondria_shifted_reference_fasta_general_indexes = [
            "${mitochondria_shifted_reference_fasta}.fai",
            sub(mitochondria_shifted_reference_fasta, ".fa(sta)?$", ".dict")
        ]

        Array[File] individual_ids
        Array[File] individual_baseline_gvcfs
        Array[File] individual_baseline_gvcf_indexes
        Array[File] individual_mitochondria_shifted_gvcfs
        Array[File] individual_mitochondria_shifted_gvcf_indexes

        String region_contig
        Int region_start
        Int region_end

        String output_prefix

        String? accelaration

        Int genomicsdb_import_batch_size = 256

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"
        String sentieon_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-sentieon:bcftools_1.11-sentieon_202010.02"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1
        String genotype_gvcfs_sentieon_docker_image = sentieon_docker_image
        Int genotype_gvcfs_sentieon_threads = 4
        Float genotype_gvcfs_sentieon_memory_gb = 16
        String vcfconcat_sentieon_docker_image = sentieon_docker_image
        Int vcfconcat_sentieon_threads = 4
        Float vcfconcat_sentieon_memory_gb = 16
        String genotype_gvcfs_cpu_docker_image = gatk_docker_image
        Int genotype_gvcfs_cpu_threads = 4
        Float genotype_gvcfs_cpu_memory_gb = 16
        String vcfconcat_cpu_docker_image = gatk_docker_image
        Int vcfconcat_cpu_threads = 4
        Float vcfconcat_cpu_memory_gb = 16

        String extract_snv_docker_image = bcftools_docker_image
        Int extract_snv_threads = 1
        Float extract_snv_memory_gb = 2
        String shift_back_docker_image = python_docker_image
        Int shift_back_threads = 2
        Float shift_back_memory_gb = 4
        String sort_vcf_docker_image = bcftools_docker_image
        Int sort_vcf_threads = 2
        Float sort_vcf_memory_gb = 4
        String compare_vcf_docker_image = bcftools_docker_image
        Int compare_vcf_threads = 4
        Float compare_vcf_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # gvcf2vcf
    # --------------------------------------------------------------------------------

    call gvcf2vcf_multi_auto.gvcf2vcf_multi_auto as step0001_baseline_gvcf2vcf { input:
        reference_fasta = baseline_reference_fasta,
        reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
        sample_ids = individual_ids,
        sample_gvcfs = individual_baseline_gvcfs,
        sample_gvcf_indexes = individual_baseline_gvcf_indexes,
        region_contig = region_contig,
        region_start = region_start,
        region_end = region_end,
        chunk_size = 1000000,
        output_prefix = "${output_prefix}.baseline",
        accelaration = accelaration,
        genomicsdb_import_batch_size = genomicsdb_import_batch_size,
        bcftools_docker_image = bcftools_docker_image,
        python_docker_image = python_docker_image,
        sentieon_docker_image = sentieon_docker_image,
        split_region_docker_image = split_region_docker_image,
        split_region_memory_gb = split_region_memory_gb,
        genotype_gvcfs_sentieon_docker_image = genotype_gvcfs_sentieon_docker_image,
        genotype_gvcfs_sentieon_threads = genotype_gvcfs_sentieon_threads,
        genotype_gvcfs_sentieon_memory_gb = genotype_gvcfs_sentieon_memory_gb,
        vcfconcat_sentieon_docker_image = vcfconcat_sentieon_docker_image,
        vcfconcat_sentieon_threads = vcfconcat_sentieon_threads,
        vcfconcat_sentieon_memory_gb = vcfconcat_sentieon_memory_gb,
        genotype_gvcfs_cpu_docker_image = genotype_gvcfs_cpu_docker_image,
        genotype_gvcfs_cpu_threads = genotype_gvcfs_cpu_threads,
        genotype_gvcfs_cpu_memory_gb = genotype_gvcfs_cpu_memory_gb,
        vcfconcat_cpu_docker_image = vcfconcat_cpu_docker_image,
        vcfconcat_cpu_threads = vcfconcat_cpu_threads,
        vcfconcat_cpu_memory_gb = vcfconcat_cpu_memory_gb
    }

    call gvcf2vcf_multi_auto.gvcf2vcf_multi_auto as step0002_shifted_gvcf2vcf { input:
        reference_fasta = mitochondria_shifted_reference_fasta,
        reference_fasta_general_indexes = mitochondria_shifted_reference_fasta_general_indexes,
        sample_ids = individual_ids,
        sample_gvcfs = individual_mitochondria_shifted_gvcfs,
        sample_gvcf_indexes = individual_mitochondria_shifted_gvcf_indexes,
        region_contig = region_contig,
        region_start = region_start,
        region_end = region_end,
        chunk_size = 1000000,
        output_prefix = "${output_prefix}.mitochondria_shifted",
        accelaration = accelaration,
        genomicsdb_import_batch_size = genomicsdb_import_batch_size,
        bcftools_docker_image = bcftools_docker_image,
        python_docker_image = python_docker_image,
        sentieon_docker_image = sentieon_docker_image,
        split_region_docker_image = split_region_docker_image,
        split_region_memory_gb = split_region_memory_gb,
        genotype_gvcfs_sentieon_docker_image = genotype_gvcfs_sentieon_docker_image,
        genotype_gvcfs_sentieon_threads = genotype_gvcfs_sentieon_threads,
        genotype_gvcfs_sentieon_memory_gb = genotype_gvcfs_sentieon_memory_gb,
        vcfconcat_sentieon_docker_image = vcfconcat_sentieon_docker_image,
        vcfconcat_sentieon_threads = vcfconcat_sentieon_threads,
        vcfconcat_sentieon_memory_gb = vcfconcat_sentieon_memory_gb,
        genotype_gvcfs_cpu_docker_image = genotype_gvcfs_cpu_docker_image,
        genotype_gvcfs_cpu_threads = genotype_gvcfs_cpu_threads,
        genotype_gvcfs_cpu_memory_gb = genotype_gvcfs_cpu_memory_gb,
        vcfconcat_cpu_docker_image = vcfconcat_cpu_docker_image,
        vcfconcat_cpu_threads = vcfconcat_cpu_threads,
        vcfconcat_cpu_memory_gb = vcfconcat_cpu_memory_gb
    }

    # --------------------------------------------------------------------------------
    # merge
    # --------------------------------------------------------------------------------

    call mitochondria_merge.mitochondria_merge as step1001_mitochondria_merge { input:
        baseline_reference_fasta = baseline_reference_fasta,
        baseline_reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
        mitochondria_shifted_reference_fasta = mitochondria_shifted_reference_fasta,
        mitochondria_shifted_reference_fasta_general_indexes = mitochondria_shifted_reference_fasta_general_indexes,
        baseline_vcf = step0001_baseline_gvcf2vcf.vcf_gz,
        baseline_vcf_index = step0001_baseline_gvcf2vcf.vcf_gz_index,
        shifted_vcf = step0002_shifted_gvcf2vcf.vcf_gz,
        shifted_vcf_index = step0002_shifted_gvcf2vcf.vcf_gz_index,
        output_prefix = output_prefix,
        bcftools_docker_image = bcftools_docker_image,
        gatk_docker_image = gatk_docker_image,
        extract_snv_docker_image = extract_snv_docker_image,
        extract_snv_threads = extract_snv_threads,
        extract_snv_memory_gb = extract_snv_memory_gb,
        shift_back_docker_image = shift_back_docker_image,
        shift_back_threads = shift_back_threads,
        shift_back_memory_gb = shift_back_memory_gb,
        sort_vcf_docker_image = sort_vcf_docker_image,
        sort_vcf_threads = sort_vcf_threads,
        sort_vcf_memory_gb = sort_vcf_memory_gb,
        compare_vcf_docker_image = compare_vcf_docker_image,
        compare_vcf_threads = compare_vcf_threads,
        compare_vcf_memory_gb = compare_vcf_memory_gb
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File merged_vcf_gz = step1001_mitochondria_merge.merged_vcf_gz
        File merged_vcf_gz_index = step1001_mitochondria_merge.merged_vcf_gz_index
        File common_baseline_vcf_gz = step1001_mitochondria_merge.common_baseline_vcf_gz
        File common_baseline_vcf_gz_index = step1001_mitochondria_merge.common_baseline_vcf_gz_index
        File common_shifted_back_vcf_gz = step1001_mitochondria_merge.common_shifted_back_vcf_gz
        File common_shifted_back_vcf_gz_index = step1001_mitochondria_merge.common_shifted_back_vcf_gz_index
        File unique_baseline_vcf_gz = step1001_mitochondria_merge.unique_baseline_vcf_gz
        File unique_baseline_vcf_gz_index = step1001_mitochondria_merge.unique_baseline_vcf_gz_index
        File unique_shifted_back_vcf_gz = step1001_mitochondria_merge.unique_shifted_back_vcf_gz
        File unique_shifted_back_vcf_gz_index = step1001_mitochondria_merge.unique_shifted_back_vcf_gz_index
    }

}
