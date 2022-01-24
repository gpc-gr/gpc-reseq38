#
#
#

version 1.0

import "./cram2gvcf_cpu.wdl"


workflow cram2gvcf_auto {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        Array[Object] regions

        String sample_id
        Int sample_ploidy
        File cram
        File cram_index
        File? bqsr_table
        Float? snv_heterozygosity
        Float? indel_heterozygosity

        String gvcf_gz__name

        String? accelaration

        String apply_bqsr_cpu_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int apply_bqsr_cpu_threads = 1
        Float apply_bqsr_cpu_memory_gb = 8
        String haplotype_caller_cpu_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int haplotype_caller_cpu_threads = 4
        Float haplotype_caller_cpu_memory_gb = 4
        String concat_gvcf_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        Int concat_gvcf_threads = 4
        Float concat_gvcf_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # cram2gvcf
    # --------------------------------------------------------------------------------

    if (!defined(accelaration)) {

        call cram2gvcf_cpu.cram2gvcf_cpu as cram2gvcf_cpu { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            regions = regions,
            sample_id = sample_id,
            sample_ploidy = sample_ploidy,
            cram = cram,
            cram_index = cram_index,
            bqsr_table = bqsr_table,
            gvcf_gz__name = gvcf_gz__name,
            snv_heterozygosity = snv_heterozygosity,
            indel_heterozygosity = indel_heterozygosity,
            apply_bqsr_cpu_docker_image = apply_bqsr_cpu_docker_image,
            apply_bqsr_cpu_threads = apply_bqsr_cpu_threads,
            apply_bqsr_cpu_memory_gb = apply_bqsr_cpu_memory_gb,
            haplotype_caller_cpu_docker_image = haplotype_caller_cpu_docker_image,
            haplotype_caller_cpu_threads = haplotype_caller_cpu_threads,
            haplotype_caller_cpu_memory_gb = haplotype_caller_cpu_memory_gb,
            concat_gvcf_docker_image = concat_gvcf_docker_image,
            concat_gvcf_threads = concat_gvcf_threads,
            concat_gvcf_memory_gb = concat_gvcf_memory_gb
        }

    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File gvcf_gz = select_first([
            cram2gvcf_cpu.gvcf_gz
        ])
        File gvcf_gz_index = select_first([
            cram2gvcf_cpu.gvcf_gz_index
        ])
    }

}
