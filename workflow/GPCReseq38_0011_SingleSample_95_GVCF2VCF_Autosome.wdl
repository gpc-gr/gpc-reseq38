#
#
#

version 1.0

import "./modules/gvcf2vcf_direct.wdl"
import "./modules/md5sum.wdl"


workflow GPCReseq38_0011_SingleSample_95_GVCF2VCF_Autosome {

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
        GVCF sample_gvcf

        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String genotype_gvcfs_docker_image = gatk_docker_image
        Int genotype_gvcfs_threads = 2
        Float genotype_gvcfs_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # gvcf2vcf
    # --------------------------------------------------------------------------------

    call gvcf2vcf_direct.gvcf2vcf_direct as step0001_gvcf2vcf_direct { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        gvcf = sample_gvcf.gvcf,
        gvcf_index = sample_gvcf.gvcf_index,
        vcf_gz__name = sub(basename(sample_gvcf.gvcf), ".g.vcf.gz$", ".vcf.gz"),
        gatk_docker_image = gatk_docker_image,
        python_docker_image = python_docker_image,
        genotype_gvcfs_docker_image = genotype_gvcfs_docker_image,
        genotype_gvcfs_threads = genotype_gvcfs_threads,
        genotype_gvcfs_memory_gb = genotype_gvcfs_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = [
            step0001_gvcf2vcf_direct.vcf_gz,
            step0001_gvcf2vcf_direct.vcf_gz_index,
        ],
        md5sum_txt__name = "${sample_id}.GPCReseq_0011_SingleSample_95_GVCF2VCF_Autosome.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File sample_vcf_gz = step0001_gvcf2vcf_direct.vcf_gz
        File sample_vcf_gz_index = step0001_gvcf2vcf_direct.vcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct GVCF {
    File gvcf
    File gvcf_index
}
