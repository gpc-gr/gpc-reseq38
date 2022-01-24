#
#
#

version 1.0


workflow fastq2cram_cpu {

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

        Boolean apply_rmdup = true

        String align_cpu_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-bwa:bwa_0.7.15-samtools_1.11"
        Int align_cpu_total_threads = 8
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
    # alignment
    # --------------------------------------------------------------------------------

    scatter (read_pair in read_pairs) {

        call step0001_align_and_sort { input:
            reference_fasta = reference_fasta,
            reference_fasta_bwa2_indexes = reference_fasta_bwa2_indexes,
            rg = read_pair.left.right,
            read1 = read_pair.right.left,
            read2 = read_pair.right.right,
            bam__name = "${read_pair.left.left}.bam",
            docker_image = align_cpu_docker_image,
            total_threads = align_cpu_total_threads,
            bwa_threads = align_cpu_bwa_mem_threads,
            samtools_threads = align_cpu_samtools_sort_threads,
            memory_gb = align_cpu_memory_gb
        }

    }

    # --------------------------------------------------------------------------------
    # rmdup & convert
    # --------------------------------------------------------------------------------

    call step0002_rmdup { input:
        source_bams = step0001_align_and_sort.bam,
        bam__name = sub(cram__name, ".cram$", ".bam"),
        rmdup_metrics__name = "${cram__name}.picard.rmdup_metrics",
        docker_image = rmdup_cpu_docker_image,
        threads = rmdup_cpu_threads,
        memory_gb = rmdup_cpu_memory_gb
    }

    call step1001_convert_bam_to_cram { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        bam = step0002_rmdup.bam,
        cram__name = cram__name,
        docker_image = bam2cram_docker_image,
        threads = bam2cram_threads,
        memory_gb = bam2cram_memory_gb
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File cram = step1001_convert_bam_to_cram.cram
        File cram_index = step1001_convert_bam_to_cram.cram_index
        File rmdup_metrics = step0002_rmdup.rmdup_metrics
    }

}


task step0001_align_and_sort {

    input {
        File reference_fasta
        Array[File] reference_fasta_bwa2_indexes
        String rg
        File read1
        File read2
        String bam__name

        String docker_image
        Int total_threads
        Int bwa_threads
        Int samtools_threads
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: total_threads
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        bwa-mem2 mem \
            -t ~{bwa_threads} \
            -R '~{rg}' \
            -T 0 \
            -Y \
            -K 10000000 \
            ~{reference_fasta} \
            ~{read1} \
            ~{read2} \
        | samtools sort \
            --threads ~{samtools_threads} \
            -l 1 \
            -o ~{bam__name} \
            -T ~{bam__name}.temp \
            -
    >>>

    output {
        File bam = "${bam__name}"
    }

}


task step0002_rmdup {

    input {
        Array[File] source_bams
        String bam__name
        String rmdup_metrics__name

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
        set -eux
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        gatk MarkDuplicates \
            ~{sep=" " prefix("--INPUT=", source_bams)} \
            --ASSUME_SORTED=true \
            --OUTPUT=~{bam__name} \
            --COMPRESSION_LEVEL=5 \
            --CREATE_INDEX=true \
            --METRICS_FILE=~{rmdup_metrics__name} \
            --VALIDATION_STRINGENCY=LENIENT \
            --USE_JDK_DEFLATER=true \
            --USE_JDK_INFLATER=true

        mv ~{sub(bam__name, ".bam$", ".bai")} ~{bam__name}.bai
    >>>

    output {
        File bam = "${bam__name}"
        File bam_index = "${bam__name}.bai"
        File rmdup_metrics = "${rmdup_metrics__name}"
    }

}


task step1001_convert_bam_to_cram {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        File bam
        String cram__name

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
        samtools view \
            --threads 8 \
            -C \
            --reference ~{reference_fasta} \
            -o ~{cram__name} \
            --write-index \
            ~{bam}
    >>>

    output {
        File cram = "${cram__name}"
        File cram_index = "${cram__name}.crai"
    }

}
