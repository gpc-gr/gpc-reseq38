#
#
#

version 1.0

import "./modules/md5sum.wdl"
import "./modules/vcf2plink.wdl"


workflow GPCReseq38_0013_Joint_18_PedigreeInference {

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

        Array[VCF] vcfs
        File sample_idtable
        String sample_sex_column_key

        Int chunk_size = 3 * 1000 * 1000    # 3Mb
        Float hwe_ge = 0.0001
        Float maf_ge = 0.01
        Float geno_le = 0.01

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String king_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-king:king_2.2.7"
        String plink_docker_image = "quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"
        String python_docker_image = "python:3.8.6-slim-buster"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1
        String make_slim_vcf_docker_image = bcftools_docker_image
        Int make_slim_vcf_threads = 4
        Float make_slim_vcf_memory_gb = 16
        String concat_vcfs_docker_image = bcftools_docker_image
        Int concat_vcfs_threads = 8
        Float concat_vcfs_memory_gb = 16
        String vcf2plink_docker_image = plink_docker_image
        Int vcf2plink_threads = 1
        Float vcf2plink_memory_gb = 90
        String update_sex_in_fam_docker_image = python_docker_image
        Float update_sex_in_fam_memory_gb = 4
        String king_inference_docker_image = king_docker_image
        Int king_inference_threads = 16
        Float king_inference_memory_gb = 90
        String update_parents_in_fam_docker_image = plink_docker_image
        Float update_parents_in_fam_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # vcf2plink
    # --------------------------------------------------------------------------------

    scatter (entry in vcfs) {
        File source_vcfs = entry.vcf
        File source_vcf_indexes = entry.vcf_index
        String contigs = entry.contig
        Int starts = entry.start
        Int ends = entry.end
    }

    call vcf2plink.vcf2plink as step0001_vcf2plink_pass_biallelic_snv { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        source_vcfs = source_vcfs,
        source_vcf_indexes = source_vcf_indexes,
        contigs = contigs,
        starts = starts,
        ends = ends,
        output_prefix = "${analysis_id}.QCed",
        pass_only = true,
        snv_only = true,
        biallelic_only = true,
        hwe_ge = hwe_ge,
        maf_ge = maf_ge,
        geno_le = geno_le,
        normalize_position_and_alleles = true,
        chunk_size = chunk_size,
        sample_idtable = sample_idtable,
        sample_sex_column_key = sample_sex_column_key,
        split_region_docker_image = split_region_docker_image,
        split_region_memory_gb = split_region_memory_gb,
        make_slim_vcf_docker_image = make_slim_vcf_docker_image,
        make_slim_vcf_threads = make_slim_vcf_threads,
        make_slim_vcf_memory_gb = make_slim_vcf_memory_gb,
        concat_vcfs_docker_image = concat_vcfs_docker_image,
        concat_vcfs_threads = concat_vcfs_threads,
        concat_vcfs_memory_gb = concat_vcfs_memory_gb,
        vcf2plink_docker_image = vcf2plink_docker_image,
        vcf2plink_threads = vcf2plink_threads,
        vcf2plink_memory_gb = vcf2plink_memory_gb,
        update_sex_in_fam_docker_image = update_sex_in_fam_docker_image,
        update_sex_in_fam_memory_gb = update_sex_in_fam_memory_gb
    }

    # --------------------------------------------------------------------------------
    # pedigree estimation
    # --------------------------------------------------------------------------------

    call step1001_king { input:
        bed = step0001_vcf2plink_pass_biallelic_snv.bed,
        bim = step0001_vcf2plink_pass_biallelic_snv.bim,
        fam = step0001_vcf2plink_pass_biallelic_snv.fam,
        output_prefix = "${analysis_id}.QCed.king",
        docker_image = king_inference_docker_image,
        threads = king_inference_threads,
        memory_gb = king_inference_memory_gb
    }

    call step1002_update_fam { input:
        original_bed = step0001_vcf2plink_pass_biallelic_snv.bed,
        original_bim = step0001_vcf2plink_pass_biallelic_snv.bim,
        original_fam = step0001_vcf2plink_pass_biallelic_snv.fam,
        king_update_ids_txt = step1001_king.king_update_ids_txt,
        king_update_parents_txt = step1001_king.king_update_parents_txt,
        docker_image = update_parents_in_fam_docker_image,
        memory_gb = update_parents_in_fam_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    Array[File] pass_snv_plink_files_selected = select_all([
        step0001_vcf2plink_pass_biallelic_snv.bed,
        step0001_vcf2plink_pass_biallelic_snv.bim,
        step1002_update_fam.result_fam
    ])

    Array[File] pedigree_inference_files_selected = select_all([
        step1001_king.king_kin,
        step1001_king.king_kin0,
        step1001_king.king_allsegs_txt,
        step1001_king.king_build_log,
        step1001_king.king_unrelated_txt,
        step1001_king.king_unrelated_to_be_removed_txt,
        step1001_king.king_update_ids_txt,
        step1001_king.king_update_parents_txt
    ])

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            pass_snv_plink_files_selected,
            pedigree_inference_files_selected
        ]),
        md5sum_txt__name = "${analysis_id}.GPCReseq38_0013_Joint_18_PedigreeInference.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        Array[File] pass_snv_plink_files = pass_snv_plink_files_selected
        Array[File] pedigree_inference_files = pedigree_inference_files_selected

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct VCF {
    File vcf
    File vcf_index
    String contig
    Int start
    Int end
}


task step1001_king {

    input {
        File bed
        File bim
        File fam
        String output_prefix

        String docker_image
        Int threads
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        king \
            -b ~{bed} \
            --bim ~{bim} \
            --fam ~{fam} \
            --kinship \
            --build \
            --unrelated \
            --degree 2 \
            --prefix ~{output_prefix}. \
            --cpus ~{threads}

        mv ~{output_prefix}..kin ~{output_prefix}.kin
        mv ~{output_prefix}..kin0 ~{output_prefix}.kin0
    >>>

    output {
        File king_kin = "${output_prefix}.kin"
        File king_kin0 = "${output_prefix}.kin0"
        File king_allsegs_txt = "${output_prefix}.allsegs.txt"
        File king_build_log = "${output_prefix}.build.log"
        File king_unrelated_txt = "${output_prefix}.unrelated.txt"
        File king_unrelated_to_be_removed_txt = "${output_prefix}.unrelated_toberemoved.txt"
        File king_update_ids_txt = "${output_prefix}.updateids.txt"
        File king_update_parents_txt = "${output_prefix}.updateparents.txt"
    }

}


task step1002_update_fam {

    input {
        File original_bed
        File original_bim
        File original_fam
        File king_update_ids_txt
        File king_update_parents_txt

        String docker_image
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        plink \
            --bed ~{original_bed} \
            --bim ~{original_bim} \
            --fam ~{original_fam} \
            --update-ids ~{king_update_ids_txt} \
            --make-just-fam \
            --out temp \
            --memory ~{floor(memory_gb * 1000 * 8/10)}

        plink \
            --bed ~{original_bed} \
            --bim ~{original_bim} \
            --fam temp.fam \
            --update-parents ~{king_update_parents_txt} \
            --make-just-fam \
            --out ~{sub(basename(original_fam), ".fam$", "")} \
            --memory ~{floor(memory_gb * 1000 * 8/10)}

        rm temp.*
    >>>

    output {
        File result_fam = "${basename(original_fam)}"
    }

}
