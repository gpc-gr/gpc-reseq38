#
#
#

version 1.0

import "./modules/md5sum.wdl"
import "./modules/mitochondria_joint.wdl"


workflow GPCReseq38_0013_Joint_02_JointGenotyping_Mitochondria {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        String analysis_id

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

        Array[GVCF] sample_gvcfs
        File target_region_list

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
    # gvcf2vcf & merge
    # --------------------------------------------------------------------------------

    scatter (entry in sample_gvcfs) {
        String individual_ids = entry.id
        File individual_baseline_gvcfs = entry.baseline_gvcf
        File individual_baseline_gvcf_indexes = entry.baseline_gvcf_index
        File individual_mitochondria_shifted_gvcfs = entry.mitochondria_shifted_gvcf
        File individual_mitochondria_shifted_gvcf_indexes = entry.mitochondria_shifted_gvcf_index
    }

    scatter (region in read_objects(target_region_list)) {

        call mitochondria_joint.mitochondria_joint as step0001_mitochondria_joint { input:
            baseline_reference_fasta = baseline_reference_fasta,
            baseline_reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
            mitochondria_shifted_reference_fasta = mitochondria_shifted_reference_fasta,
            mitochondria_shifted_reference_fasta_general_indexes = mitochondria_shifted_reference_fasta_general_indexes,
            individual_ids = individual_ids,
            individual_baseline_gvcfs = individual_baseline_gvcfs,
            individual_baseline_gvcf_indexes = individual_baseline_gvcf_indexes,
            individual_mitochondria_shifted_gvcfs = individual_mitochondria_shifted_gvcfs,
            individual_mitochondria_shifted_gvcf_indexes = individual_mitochondria_shifted_gvcf_indexes,
            region_contig = region.contig,
            region_start = region.start,
            region_end = region.end,
            output_prefix = "${analysis_id}.${region.name}",
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

    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            step0001_mitochondria_joint.merged_vcf_gz,
            step0001_mitochondria_joint.merged_vcf_gz_index,
            step0001_mitochondria_joint.common_baseline_vcf_gz,
            step0001_mitochondria_joint.common_baseline_vcf_gz_index,
            step0001_mitochondria_joint.common_shifted_back_vcf_gz,
            step0001_mitochondria_joint.common_shifted_back_vcf_gz_index,
            step0001_mitochondria_joint.unique_baseline_vcf_gz,
            step0001_mitochondria_joint.unique_baseline_vcf_gz_index,
            step0001_mitochondria_joint.unique_shifted_back_vcf_gz,
            step0001_mitochondria_joint.unique_shifted_back_vcf_gz_index
        ]),
        md5sum_txt__name = "${analysis_id}.GPCReseq38_0013_Joint_02_JointGenotyping_Mitochondria.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        Array[File] merged_vcf_gz = step0001_mitochondria_joint.merged_vcf_gz
        Array[File] merged_vcf_gz_index = step0001_mitochondria_joint.merged_vcf_gz_index
        Array[File] common_baseline_vcf_gz = step0001_mitochondria_joint.common_baseline_vcf_gz
        Array[File] common_baseline_vcf_gz_index = step0001_mitochondria_joint.common_baseline_vcf_gz_index
        Array[File] common_shifted_back_vcf_gz = step0001_mitochondria_joint.common_shifted_back_vcf_gz
        Array[File] common_shifted_back_vcf_gz_index = step0001_mitochondria_joint.common_shifted_back_vcf_gz_index
        Array[File] unique_baseline_vcf_gz = step0001_mitochondria_joint.unique_baseline_vcf_gz
        Array[File] unique_baseline_vcf_gz_index = step0001_mitochondria_joint.unique_baseline_vcf_gz_index
        Array[File] unique_shifted_back_vcf_gz = step0001_mitochondria_joint.unique_shifted_back_vcf_gz
        Array[File] unique_shifted_back_vcf_gz_index = step0001_mitochondria_joint.unique_shifted_back_vcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct GVCF {
    String id
    File baseline_gvcf
    File baseline_gvcf_index
    File mitochondria_shifted_gvcf
    File mitochondria_shifted_gvcf_index
}
