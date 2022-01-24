#
#
#

version 1.0

import "./modules/gvcf2vcf_direct.wdl"
import "./modules/md5sum.wdl"
import "./modules/mitochondria_merge.wdl"


workflow GPCReseq38_0011_SingleSample_98_GVCF2VCF_Mitochondria {

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

        String sample_id
        String sample_id_with_suffix
        GVCF sample_baseline_gvcf
        GVCF sample_shifted_gvcf

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String genotype_gvcfs_docker_image = gatk_docker_image
        Int genotype_gvcfs_threads = 2
        Float genotype_gvcfs_memory_gb = 4

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

    call gvcf2vcf_direct.gvcf2vcf_direct as step0001_baseline_gvcf2vcf { input:
        reference_fasta = baseline_reference_fasta,
        reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
        gvcf = sample_baseline_gvcf.gvcf,
        gvcf_index = sample_baseline_gvcf.gvcf_index,
        vcf_gz__name = sub(basename(sample_baseline_gvcf.gvcf), ".g.vcf.gz$", ".all.vcf.gz"),
        gatk_docker_image = gatk_docker_image,
        genotype_gvcfs_docker_image = genotype_gvcfs_docker_image,
        genotype_gvcfs_threads = genotype_gvcfs_threads,
        genotype_gvcfs_memory_gb = genotype_gvcfs_memory_gb
    }

    call gvcf2vcf_direct.gvcf2vcf_direct as step0002_shifted_gvcf2vcf { input:
        reference_fasta = mitochondria_shifted_reference_fasta,
        reference_fasta_general_indexes = mitochondria_shifted_reference_fasta_general_indexes,
        gvcf = sample_shifted_gvcf.gvcf,
        gvcf_index = sample_shifted_gvcf.gvcf_index,
        vcf_gz__name = sub(basename(sample_shifted_gvcf.gvcf), ".g.vcf.gz$", ".vcf.gz"),
        gatk_docker_image = gatk_docker_image,
        genotype_gvcfs_docker_image = genotype_gvcfs_docker_image,
        genotype_gvcfs_threads = genotype_gvcfs_threads,
        genotype_gvcfs_memory_gb = genotype_gvcfs_memory_gb
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
        output_prefix = sample_id_with_suffix,
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
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = [
            step1001_mitochondria_merge.merged_vcf_gz,
            step1001_mitochondria_merge.merged_vcf_gz_index,
            step1001_mitochondria_merge.common_baseline_vcf_gz,
            step1001_mitochondria_merge.common_baseline_vcf_gz_index,
            step1001_mitochondria_merge.common_shifted_back_vcf_gz,
            step1001_mitochondria_merge.common_shifted_back_vcf_gz_index,
            step1001_mitochondria_merge.unique_baseline_vcf_gz,
            step1001_mitochondria_merge.unique_baseline_vcf_gz_index,
            step1001_mitochondria_merge.unique_shifted_back_vcf_gz,
            step1001_mitochondria_merge.unique_shifted_back_vcf_gz_index
        ],
        md5sum_txt__name = "${sample_id_with_suffix}.GPCReseq_0011_SingleSample_98_GVCF2VCF_Mitochondria.md5sum.txt",
        docker_image = python_docker_image
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

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct GVCF {
    File gvcf
    File gvcf_index
}
