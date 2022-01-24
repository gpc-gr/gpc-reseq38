#
#
#

version 1.0

import "./modules/gvcf2vcf_multi_auto.wdl"
import "./modules/md5sum.wdl"


workflow GPCReseq38_0013_Joint_01_JointGenotyping {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        String analysis_id

        File reference_fasta
        Array[File] reference_fasta_general_indexes = [
            "${reference_fasta}.fai",
            sub(reference_fasta, ".fa(sta)?$", ".dict")
        ]

        Array[GVCF] sample_gvcfs

        File target_region_list
        Int chunk_size = 3 * 1000 * 1000    # 3Mb
        Int chunk_padding = 3000            # 3kb

        Int genomicsdb_import_batch_size = 256

        String? accelaration

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
        Int genotype_gvcfs_cpu_threads = 3
        Float genotype_gvcfs_cpu_memory_gb = 6
        String vcfconcat_cpu_docker_image = bcftools_docker_image
        Int vcfconcat_cpu_threads = 4
        Float vcfconcat_cpu_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # joint genotyping
    # --------------------------------------------------------------------------------

    scatter (entry in sample_gvcfs) {

        String individual_id = entry.id
        File individual_gvcf = entry.gvcf
        File individual_gvcf_index = entry.gvcf_index

    }

    scatter (region in read_objects(target_region_list)) {

        call gvcf2vcf_multi_auto.gvcf2vcf_multi_auto as step0001_gvcf2vcf { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            sample_ids = individual_id,
            sample_gvcfs = individual_gvcf,
            sample_gvcf_indexes = individual_gvcf_index,
            region_contig = region.contig,
            region_start = region.start,
            region_end = region.end,
            chunk_size = chunk_size,
            chunk_padding = chunk_padding,
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
            step0001_gvcf2vcf.vcf_gz,
            step0001_gvcf2vcf.vcf_gz_index
        ]),
        md5sum_txt__name = "${analysis_id}.GPCReseq38_0013_Joint_01_JointGenotyping.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        Array[File] vcf_gz = step0001_gvcf2vcf.vcf_gz
        Array[File] vcf_gz_index = step0001_gvcf2vcf.vcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct GVCF {
    String id
    File gvcf
    File gvcf_index
}
