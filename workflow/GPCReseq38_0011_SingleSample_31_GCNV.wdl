#
#
#

version 1.0

import "./modules/md5sum.wdl"


workflow GPCReseq38_0011_SingleSample_31_GCNV {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes = [
            "${reference_fasta}.fai",
            sub(reference_fasta, ".fa(sta)?$", ".dict")
        ]

        File interval_list

        String sample_id
        File sample_cram
        File sample_cram_index = "${sample_cram}.crai"

        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String collect_read_counts_docker_image = gatk_docker_image
        Int collect_read_counts_threads = 1
        Float collect_read_counts_memory_gb = 16
    }

    # --------------------------------------------------------------------------------
    # main
    # --------------------------------------------------------------------------------

    call step0001_collect_read_counts { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        interval_list = interval_list,
        cram = sample_cram,
        cram_index = sample_cram_index,
        hdf5__name = "${basename(sample_cram)}.gatk.collect_read_counts.hdf5",
        docker_image = gatk_docker_image,
        threads = collect_read_counts_threads,
        memory_gb = collect_read_counts_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = [
            step0001_collect_read_counts.hdf5
        ],
        md5sum_txt__name = "${sample_id}.GPCReseq_0011_SingleSample_31_GCNV.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File read_count_hdf5 = step0001_collect_read_counts.hdf5
        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


task step0001_collect_read_counts {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        File interval_list
        File cram
        File cram_index
        String hdf5__name

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
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        gatk CollectReadCounts \
            --reference ~{reference_fasta} \
            --intervals ~{interval_list} \
            --interval-merging-rule OVERLAPPING_ONLY \
            --input ~{cram} \
            --format HDF5 \
            --output ~{hdf5__name}
    >>>

    output {
        File hdf5 = "${hdf5__name}"
    }

}
