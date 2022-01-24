#
#
#

version 1.0

import "./modules/cram2gvcf_auto.wdl"
import "./modules/md5sum.wdl"


workflow GPCReseq38_0011_SingleSample_12_chrXY_PAR2 {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes = [
            "${reference_fasta}.fai",
            sub(reference_fasta, ".fa(sta)?$", ".dict")
        ]

        String sample_id
        File sample_cram
        File sample_cram_index = "${sample_cram}.crai"
        File? sample_bqsr_table

        String? accelaration

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String apply_bqsr_cpu_docker_image = gatk_docker_image
        Int apply_bqsr_cpu_threads = 1
        Float apply_bqsr_cpu_memory_gb = 8
        String haplotype_caller_cpu_docker_image = gatk_docker_image
        Int haplotype_caller_cpu_threads = 4
        Float haplotype_caller_cpu_memory_gb = 4
        String concat_gvcf_docker_image = bcftools_docker_image
        Int concat_gvcf_threads = 4
        Float concat_gvcf_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # variant call
    # --------------------------------------------------------------------------------

    String sample_id_with_suffix = sample_id + (if (defined(sample_bqsr_table)) then ".BQSR" else ".noBQSR")
    Array[VariantCallConfig] configs = [
        {
            "region_list": "${reference_fasta}.regions.chrXY_PAR2.chrX_nonPAR.tsv",
            "ploidy": 1,
            "gvcf_gz__name": "${sample_id_with_suffix}.chrXY_PAR2.chrX_nonPAR.male.g.vcf.gz"
        },
        {
            "region_list": "${reference_fasta}.regions.chrXY_PAR2.chrX_nonPAR.tsv",
            "ploidy": 2,
            "gvcf_gz__name": "${sample_id_with_suffix}.chrXY_PAR2.chrX_nonPAR.female.g.vcf.gz"
        },
        {
            "region_list": "${reference_fasta}.regions.chrXY_PAR2.chrX_PAR.tsv",
            "ploidy": 2,
            "gvcf_gz__name": "${sample_id_with_suffix}.chrXY_PAR2.chrX_PAR.g.vcf.gz"
        },
        {
            "region_list": "${reference_fasta}.regions.chrXY_PAR2.chrY_nonPAR.tsv",
            "ploidy": 1,
            "gvcf_gz__name": "${sample_id_with_suffix}.chrXY_PAR2.chrY_nonPAR.male.g.vcf.gz"
        }
    ]

    scatter (config in configs) {

        call cram2gvcf_auto.cram2gvcf_auto as step0001_cram2gvcf { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            regions = read_objects(config.region_list),
            sample_id = sample_id,
            sample_ploidy = config.ploidy,
            cram = sample_cram,
            cram_index = sample_cram_index,
            bqsr_table = sample_bqsr_table,
            gvcf_gz__name = config.gvcf_gz__name,
            accelaration = accelaration,
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
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            step0001_cram2gvcf.gvcf_gz,
            step0001_cram2gvcf.gvcf_gz_index
        ]),
        md5sum_txt__name = "${sample_id_with_suffix}.GPCReseq_0011_SingleSample_12_chrXY_PAR2.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        Array[File] gvcf_gz = step0001_cram2gvcf.gvcf_gz
        Array[File] gvcf_gz_index = step0001_cram2gvcf.gvcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct VariantCallConfig {
    File region_list
    Int ploidy
    String gvcf_gz__name
}
