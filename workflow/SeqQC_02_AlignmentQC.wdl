#
#
#

version 1.0

import "./modules/multiqc.wdl"


workflow SeqQC_02_AlignmentQC {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes = [
            "${reference_fasta}.fai",
            sub(reference_fasta, ".fa(sta)?$", ".dict")
        ]
        Array[File] reference_fasta_bwa2_indexes = [
            "${reference_fasta}.0123",
            "${reference_fasta}.amb",
            "${reference_fasta}.ann",
            "${reference_fasta}.bwt.2bit.64",
            "${reference_fasta}.pac",
            "${reference_fasta}.alt"
        ]

        String flowcell
        Array[ReadPair] read_pairs

        Int? num_reads_to_align = 2000000

        String align_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-bwamem2:bwamem2_2.2.1-samtools_1.11"
        Int align_bwa_mem_threads = 8
        Int align_samtools_sort_threads = 4
        Float align_memory_gb = 32

        String picard_insert_size_metrics_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int picard_insert_size_metrics_memory_gb = 8

        String picard_alignment_summary_metrics_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int picard_alignment_summary_metrics_memory_gb = 8

        String multiqc_docker_image = "quay.io/biocontainers/multiqc:1.9--py_1"
        Int multiqc_threads = 2
        Float multiqc_memory_gb = 48
    }

    # --------------------------------------------------------------------------------
    # align & collect metrics
    # --------------------------------------------------------------------------------

    scatter (read_pair in read_pairs) {

        call step0001_align { input:
            reference_fasta = reference_fasta,
            reference_fasta_bwa2_indexes = reference_fasta_bwa2_indexes,
            read1 = read_pair.read1,
            read2 = read_pair.read2,
            num_reads_to_align = num_reads_to_align,
            docker_image = align_docker_image,
            bwa_mem_threads = align_bwa_mem_threads,
            samtools_sort_threads = align_samtools_sort_threads,
            memory_gb = align_memory_gb
        }

        call step0002_picard_insert_size_metrics { input:
            bam = step0001_align.bam,
            picard_insert_size_metrics__name = "${read_pair.id}.picard.insert_size_metrics",
            docker_image = picard_insert_size_metrics_docker_image,
            memory_gb = picard_insert_size_metrics_memory_gb
        }

        if (!defined(num_reads_to_align)) {

            call step0003_picard_alignment_summary_metrics { input:
                bam = step0001_align.bam,
                picard_alignment_summary_metrics__name = "${read_pair.id}.picard.alignment_summary_metrics",
                docker_image = picard_alignment_summary_metrics_docker_image,
                memory_gb = picard_alignment_summary_metrics_memory_gb
            }

        }

    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    Array[File] picard_metrics_files = flatten([
        step0002_picard_insert_size_metrics.picard_insert_size_metrics,
        select_all(step0003_picard_alignment_summary_metrics.picard_alignment_summary_metrics)
    ])

    call multiqc.report as step1001_multiqc { input:
        sources = picard_metrics_files,
        output_prefix = "${flowcell}.02_AlignmentQC.multiqc",
        generate_plots_zip = true,
        threads = multiqc_threads,
        memory_gb = multiqc_memory_gb
    }

    call step1002_combine_picard_metrics { input:
        sources = picard_metrics_files,
        zip__name = "${flowcell}.02_AlignmentQC.picard.metrics.zip",
        docker_image = multiqc_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File multiqc_html = step1001_multiqc.html
        File multiqc_zip = step1001_multiqc.zip
        File multiqc_plots_zip = step1001_multiqc.plots_zip

        File picard_metrics_zip = step1001_multiqc.zip
    }

}


struct ReadPair {
    String id
    File read1
    File read2
}


task step0001_align {

    input {
        File reference_fasta
        Array[File] reference_fasta_bwa2_indexes
        File read1
        File read2
        Int? num_reads_to_align
        String bam__name

        String docker_image
        Int bwa_mem_threads
        Int samtools_sort_threads
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: bwa_mem_threads + samtools_sort_threads
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        #
        if [ ~{default="-1" num_reads_to_align} -gt 0 ]; then
            gzip -dc ~{read1} | head -n ~{num_reads_to_align * 4} | gzip -c > ~{bam__name}_R1.fastq.gz || true
            if [ $(gzip -dc ~{bam__name}_R1.fastq.gz | wc -l) -eq 0 ]; then
                exit 1
            fi

            gzip -dc ~{read2} | head -n ~{num_reads_to_align * 4} | gzip -c > ~{bam__name}_R2.fastq.gz || true
            if [ $(gzip -dc ~{bam__name}_R2.fastq.gz | wc -l) -eq 0 ]; then
                exit 1
            fi

        else
            ln -s ~{read1} ~{bam__name}_R1.fastq.gz
            ln -s ~{read2} ~{bam__name}_R2.fastq.gz
        fi

        #
        bwa-mem2 mem \
            -t ~{bwa_mem_threads} \
            -T 0 \
            -Y \
            -K 10000000 \
            ~{reference_fasta} \
            ~{bam__name}_R1.fastq.gz \
            ~{bam__name}_R2.fastq.gz \
        | samtools sort \
            --threads ~{samtools_sort_threads} \
            -l 1 \
            -o ~{bam__name} \
            -T ~{bam__name}.temp \
            -

        #
        rm ~{bam__name}_R1.fastq.gz ~{bam__name}_R2.fastq.gz
    >>>

    output {
        File bam = "${bam__name}"
    }

}


task step0002_picard_insert_size_metrics {

    input {
        File bam
        File picard_insert_size_metrics__name

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
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        gatk CollectInsertSizeMetrics \
            --INPUT ~{bam} \
            --OUTPUT ~{picard_insert_size_metrics__name} \
            --Histogram_FILE /dev/null
    >>>

    output {
        File picard_insert_size_metrics = "${picard_insert_size_metrics__name}"
    }

}


task step0003_picard_alignment_summary_metrics {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        File bam
        File picard_alignment_summary_metrics__name

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
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        gatk CollectAlignmentSummaryMetrics \
            --REFERENCE_SEQUENCE ~{reference_fasta} \
            --INPUT ~{bam} \
            --OUTPUT ~{picard_alignment_summary_metrics__name}
    >>>

    output {
        File picard_alignment_summary_metrics = "${picard_alignment_summary_metrics__name}"
    }

}


task step1002_combine_picard_metrics {

    input {
        Array[File] sources
        String zip__name

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

        mkdir -p picard_metrics
        cat ~{write_lines(sources)} | xargs -n1 -I%p cp %p picard_metrics

        python -c 'import shutil; shutil.make_archive("picard_metrics", "zip", "picard_metrics")'
        mv picard_metrics.zip ~{zip__name}
    >>>

    output {
        File zip = "${zip__name}"
    }

}
