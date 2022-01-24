#
#
#

version 1.0

import "./fastq2cram_cpu.wdl"


workflow fastq2cram_auto {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        Array[File] reference_fasta_bwa2_indexes

        String sample_id
        Array[Pair[Pair[String, String], Pair[File, File]]] read_pairs
        String cram__name

        String? accelaration

        String align_cpu_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-bwamem2:bwamem2_2.2.1-samtools_1.11"
        Int align_cpu_total_threads = 10
        Int align_cpu_bwa_mem_threads = 8
        Int align_cpu_samtools_sort_threads = 4
        Float align_cpu_memory_gb = 32
        String rmdup_cpu_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int rmdup_cpu_threads = 1
        Float rmdup_cpu_memory_gb = 16

        String bam2cram_docker_image = "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        Int bam2cram_threads = 8
        Float bam2cram_memory_gb = 8
    }

    # --------------------------------------------------------------------------------
    # fastq2cram
    # --------------------------------------------------------------------------------

    if (!defined(accelaration)) {

        call fastq2cram_cpu.fastq2cram_cpu as fastq2cram_cpu  { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            reference_fasta_bwa2_indexes = reference_fasta_bwa2_indexes,
            sample_id = sample_id,
            read_pairs = read_pairs,
            cram__name = cram__name,
            align_cpu_docker_image = align_cpu_docker_image,
            align_cpu_total_threads = align_cpu_total_threads,
            align_cpu_bwa_mem_threads = align_cpu_bwa_mem_threads,
            align_cpu_samtools_sort_threads = align_cpu_samtools_sort_threads,
            align_cpu_memory_gb = align_cpu_memory_gb,
            rmdup_cpu_docker_image = rmdup_cpu_docker_image,
            rmdup_cpu_threads = rmdup_cpu_threads,
            rmdup_cpu_memory_gb = rmdup_cpu_memory_gb,
            bam2cram_docker_image = bam2cram_docker_image,
            bam2cram_threads = bam2cram_threads,
            bam2cram_memory_gb = bam2cram_memory_gb
        }

    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File cram = select_first([
            fastq2cram_cpu.cram
        ])
        File cram_index = select_first([
            fastq2cram_cpu.cram_index
        ])
        File rmdup_metrics = select_first([
            fastq2cram_cpu.rmdup_metrics
        ])
    }

}
