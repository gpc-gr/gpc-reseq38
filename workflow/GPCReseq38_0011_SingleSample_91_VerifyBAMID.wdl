#
#
#

version 1.0

import "./modules/md5sum.wdl"


workflow GPCReseq38_0011_SingleSample_91_VerifyBAMID {

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

        File verifybamid2_ud
        File verifybamid2_bed
        File verifybamid2_mean

        String python_docker_image = "python:3.8.6-slim-buster"
        String verifybamid2_image = "quay.io/biocontainers/verifybamid2:1.0.6--he56e5df_0"

        String verifybamid2_threads = 4
        String verifybamid2_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # contamination check
    # --------------------------------------------------------------------------------

    call step0001_verify_bam_id { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        verifybamid2_ud = verifybamid2_ud,
        verifybamid2_bed = verifybamid2_bed,
        verifybamid2_mean = verifybamid2_mean,
        cram = sample_cram,
        cram_index = sample_cram_index,
        out__name = "${basename(sample_cram)}.verifyBAMID2.out",
        selfsm__name = "${basename(sample_cram)}.verifyBAMID2.selfSM",
        docker_image = verifybamid2_image,
        threads = verifybamid2_threads,
        memory_gb = verifybamid2_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = [
            step0001_verify_bam_id.out,
            step0001_verify_bam_id.selfsm
        ],
        md5sum_txt__name = "${sample_id}.GPCReseq38_SingleSample_91_VerifyBAMID.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File verifybamid2_out = step0001_verify_bam_id.out
        File verifybamid2_selfsm = step0001_verify_bam_id.selfsm

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


task step0001_verify_bam_id {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        File verifybamid2_ud
        File verifybamid2_bed
        File verifybamid2_mean

        File cram
        File cram_index

        String out__name
        String selfsm__name

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

        /usr/local/share/verifybamid2-1.0.6-0/VerifyBamID \
            --Reference ~{reference_fasta} \
            --UDPath ~{verifybamid2_ud} \
            --BedPath ~{verifybamid2_bed} \
            --MeanPath ~{verifybamid2_mean} \
            --BamFile ~{cram} \
            --Seed 12345 \
            --NumThread ~{threads}

        mv result.out ~{out__name}
        mv result.selfSM ~{selfsm__name}
    >>>

    output {
        File out = "${out__name}"
        File selfsm = "${selfsm__name}"
    }

}
