#
#
#

version 1.0

import "./modules/md5sum.wdl"
import "./vendor/sv-pipeline-c49150e7/scripts/Pre_Merge_SV_per_sample.wdl"


workflow GPCReseq38_0011_SingleSample_21_SVPreMerge {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        File reference_fasta_fai = "${reference_fasta}.fai"
        File reference_fasta_cache_tar_gz = "${reference_fasta}.samtools.cache.tar.gz"

        File? call_regions_bed
        File? call_regions_bed_index
        File exclude_regions

        String sample_id
        File sample_cram
        File sample_cram_index = "${sample_cram}.crai"
        File sample_cram_suffix = ".base.cram"

        String python_docker_image = "python:3.8.6-slim-buster"
    }

    # --------------------------------------------------------------------------------
    # runs SV pipeline
    # --------------------------------------------------------------------------------

    call Pre_Merge_SV_per_sample.Pre_Merge_SV_Per_Sample as step0001_pre_merge_sv_per_sample { input:
        aligned_cram = sample_cram,
        aligned_cram_index = sample_cram_index,
        ref_fasta = reference_fasta,
        ref_fasta_index = reference_fasta_fai,
        ref_cache = reference_fasta_cache_tar_gz,
        call_regions_bed = call_regions_bed,
        call_regions_bed_index = call_regions_bed_index,
        exclude_regions = exclude_regions,
        aligned_cram_suffix = sample_cram_suffix,
        preemptible_tries = 3
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = [
            step0001_pre_merge_sv_per_sample.manta_vcf,
            step0001_pre_merge_sv_per_sample.manta_tbi,
            step0001_pre_merge_sv_per_sample.manta_original_vcf,
            step0001_pre_merge_sv_per_sample.manta_original_tbi,
            step0001_pre_merge_sv_per_sample.cnvnator_cn_hist_root,
            step0001_pre_merge_sv_per_sample.cnvnator_output_cn_txt,
            step0001_pre_merge_sv_per_sample.cnvnator_cn_bed,
            step0001_pre_merge_sv_per_sample.smoove_vcf,
            step0001_pre_merge_sv_per_sample.smoove_csi
        ],
        md5sum_txt__name = "${sample_id}.GPCReseq_0011_SingleSample_21_SVPreMerge.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File manta_vcf = step0001_pre_merge_sv_per_sample.manta_vcf
        File manta_tbi = step0001_pre_merge_sv_per_sample.manta_tbi
        File manta_original_vcf = step0001_pre_merge_sv_per_sample.manta_original_vcf
        File manta_original_tbi = step0001_pre_merge_sv_per_sample.manta_original_tbi
        File cnvnator_cn_hist_root = step0001_pre_merge_sv_per_sample.cnvnator_cn_hist_root
        File cnvnator_output_cn_txt = step0001_pre_merge_sv_per_sample.cnvnator_output_cn_txt
        File cnvnator_cn_bed = step0001_pre_merge_sv_per_sample.cnvnator_cn_bed
        File smoove_vcf = step0001_pre_merge_sv_per_sample.smoove_vcf
        File smoove_csi = step0001_pre_merge_sv_per_sample.smoove_csi

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}
