#
#
#

version 1.0

import "./modules/md5sum.wdl"


workflow GPCReseq38_0011_SingleSample_92_DepthDetail {

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

        String python_docker_image = "python:3.8.6-slim-buster"
        String mosdepth_image = "quay.io/biocontainers/mosdepth:0.3.1--h01d7912_2"

        Int mosdepth_threads = 4
        Float mosdepth_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # mosdepth
    # --------------------------------------------------------------------------------

    call mosdepth as step0001_mosdepth_mapq00 { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        window = 500,
        mapq = 0,
        cram = sample_cram,
        cram_index = sample_cram_index,
        summary_txt__name = "${basename(sample_cram)}.mosdepth.mapq00.summary.txt",
        docker_image = mosdepth_image,
        threads = mosdepth_threads,
        memory_gb = mosdepth_memory_gb
    }

    call mosdepth as step0002_mosdepth_mapq20 { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        window = 500,
        mapq = 20,
        cram = sample_cram,
        cram_index = sample_cram_index,
        summary_txt__name = "${basename(sample_cram)}.mosdepth.mapq20.summary.txt",
        docker_image = mosdepth_image,
        threads = mosdepth_threads,
        memory_gb = mosdepth_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = [
            step0001_mosdepth_mapq00.summary_txt,
            step0002_mosdepth_mapq20.summary_txt,
        ],
        md5sum_txt__name = "${sample_id}.GPCReseq38_SingleSample_92_DepthDetail.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File mosdepth_summary_txt_mapq00 = step0001_mosdepth_mapq00.summary_txt
        File mosdepth_summary_txt_mapq20 = step0002_mosdepth_mapq20.summary_txt

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


task mosdepth {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        Int window
        Int mapq

        File cram
        File cram_index

        String summary_txt__name

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

        mosdepth \
            --threads ~{threads} \
            --fasta ~{reference_fasta} \
            --no-per-base \
            --fast-mode \
            --by ~{window} \
            --mapq ~{mapq} \
            out \
            ~{cram}

        mv out.mosdepth.summary.txt ~{summary_txt__name}
        rm out.mosdepth.global.dist.txt
        rm out.mosdepth.region.dist.txt
        rm out.regions.bed.gz
        rm out.regions.bed.gz.csi
    >>>

    output {
        File summary_txt = "${summary_txt__name}"
    }

}
