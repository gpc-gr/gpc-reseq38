#
#
#

version 1.0

import "./modules/gvcf2vcf_direct.wdl"
import "./modules/md5sum.wdl"
import "./modules/vcfconcat.wdl"


workflow GPCReseq38_0011_SingleSample_97_GVCF2VCF_chrXY_PAR3 {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File chrXY_PAR3_reference_fasta
        Array[File] chrXY_PAR3_reference_fasta_general_indexes = [
            "${chrXY_PAR3_reference_fasta}.fai",
            sub(chrXY_PAR3_reference_fasta, ".fa(sta)?$", ".dict")
        ]

        String sample_id
        String sample_id_with_suffix
        Array[GVCF] sample_gvcfs

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String genotype_gvcfs_docker_image = gatk_docker_image
        Int genotype_gvcfs_threads = 2
        Float genotype_gvcfs_memory_gb = 4

        String concat_vcfs_docker_image = bcftools_docker_image
        Int concat_vcfs_threads = 4
        Float concat_vcfs_memory_gb = 1
    }

    # --------------------------------------------------------------------------------
    # gvcf2vcf
    # --------------------------------------------------------------------------------

    scatter (entry in sample_gvcfs) {

        call gvcf2vcf_direct.gvcf2vcf_direct as step0001_gvcf2vcf_direct { input:
            reference_fasta = chrXY_PAR3_reference_fasta,
            reference_fasta_general_indexes = chrXY_PAR3_reference_fasta_general_indexes,
            gvcf = entry.gvcf,
            gvcf_index = entry.gvcf_index,
            vcf_gz__name = sub(basename(entry.gvcf), ".g.vcf.gz$", ".vcf.gz"),
            gatk_docker_image = gatk_docker_image,
            python_docker_image = python_docker_image,
            genotype_gvcfs_docker_image = genotype_gvcfs_docker_image,
            genotype_gvcfs_threads = genotype_gvcfs_threads,
            genotype_gvcfs_memory_gb = genotype_gvcfs_memory_gb
        }

    }

    call vcfconcat.vcfconcat as step0002_concat_vcfs { input:
        source_vcfs = step0001_gvcf2vcf_direct.vcf_gz,
        source_vcf_indexes = step0001_gvcf2vcf_direct.vcf_gz_index,
        result_vcf_gz__name = "${sample_id_with_suffix}.chrXY_PAR3.vcf.gz",
        apply_sort = true,
        source_vcfs_have_same_header = false,
        create_index = true,
        docker_image = concat_vcfs_docker_image,
        threads = concat_vcfs_threads,
        memory_gb = concat_vcfs_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = [
            step0002_concat_vcfs.result_vcf_gz,
            step0002_concat_vcfs.result_vcf_gz_index,
        ],
        md5sum_txt__name = "${sample_id_with_suffix}.GPCReseq_0011_SingleSample_97_GVCF2VCF_chrXY_PAR3.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File sample_vcf_gz = step0002_concat_vcfs.result_vcf_gz
        File sample_vcf_gz_index = step0002_concat_vcfs.result_vcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct GVCF {
    File gvcf
    File gvcf_index
}
