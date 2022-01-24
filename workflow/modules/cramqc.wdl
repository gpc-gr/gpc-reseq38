#
#
#

version 1.0


workflow cramqc {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        Array[String] wgs_metrics_regions

        File cram
        File cram_index

        String picard_multiple_metrics_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int picard_multiple_metrics_threads = 1
        Float picard_multiple_metrics_memory_gb = 32
        String picard_wgs_metrics_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int picard_wgs_metrics_threads = 1
        Float picard_wgs_metrics_memory_gb = 32
        String samtools_idxstats_docker_image = "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        Int samtools_idxstats_threads = 8
        Float samtools_idxstats_memory_gb = 8
        String samtools_flagstat_docker_image = "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        Int samtools_flagstat_threads = 8
        Float samtools_flagstat_memory_gb = 8
    }

    # --------------------------------------------------------------------------------
    # picard
    # --------------------------------------------------------------------------------

    call step0001_picard_multiple_metrics { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        cram = cram,
        cram_index = cram_index,
        output_prefix = "${basename(cram)}.picard",
        docker_image = picard_multiple_metrics_docker_image,
        threads = picard_multiple_metrics_threads,
        memory_gb = picard_multiple_metrics_memory_gb
    }

    scatter (region in wgs_metrics_regions) {

        call step0002_picard_wgs_metrics { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            region_name = region,
            interval_list = "${reference_fasta}.regions.${region}.interval_list",
            cram = cram,
            cram_index = cram_index,
            wgs_metrics__name = "${basename(cram)}.picard.wgs_metrics.${region}",
            docker_image = picard_wgs_metrics_docker_image,
            threads = picard_wgs_metrics_threads,
            memory_gb = picard_wgs_metrics_memory_gb
        }

    }

    # --------------------------------------------------------------------------------
    # samtools
    # --------------------------------------------------------------------------------

    call step1001_samtools_idxstats { input:
        cram = cram,
        cram_index = cram_index,
        idxstats__name = "${basename(cram)}.samtools.idxstats",
        docker_image = samtools_idxstats_docker_image,
        threads = samtools_idxstats_threads,
        memory_gb = samtools_idxstats_memory_gb
    }

    call step1002_samtools_flagstat { input:
        cram = cram,
        cram_index = cram_index,
        flagstat__name = "${basename(cram)}.samtools.flagstat",
        docker_image = samtools_flagstat_docker_image,
        threads = samtools_flagstat_threads,
        memory_gb = samtools_flagstat_memory_gb
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File picard_alignment_summary_metrics = step0001_picard_multiple_metrics.alignment_summary_metrics
        File picard_insert_size_metrics = step0001_picard_multiple_metrics.insert_size_metrics
        File picard_insert_size_histogram_pdf = step0001_picard_multiple_metrics.insert_size_histogram_pdf
        File picard_quality_distribution_metrics = step0001_picard_multiple_metrics.quality_distribution_metrics
        File picard_quality_distribution_pdf = step0001_picard_multiple_metrics.quality_distribution_pdf
        File picard_quality_by_cycle_metrics = step0001_picard_multiple_metrics.quality_by_cycle_metrics
        File picard_quality_by_cycle_pdf = step0001_picard_multiple_metrics.quality_by_cycle_pdf
        File picard_base_distribution_by_cycle_metrics = step0001_picard_multiple_metrics.base_distribution_by_cycle_metrics
        File picard_base_distribution_by_cycle_pdf = step0001_picard_multiple_metrics.base_distribution_by_cycle_pdf
        File picard_gc_bias_pdf = step0001_picard_multiple_metrics.gc_bias_pdf
        File picard_gc_bias_summary_metrics = step0001_picard_multiple_metrics.gc_bias_summary_metrics
        File picard_gc_bias_detail_metrics = step0001_picard_multiple_metrics.gc_bias_detail_metrics
        File picard_bait_bias_summary_metrics = step0001_picard_multiple_metrics.bait_bias_summary_metrics
        File picard_bait_bias_detail_metrics = step0001_picard_multiple_metrics.bait_bias_detail_metrics
        File picard_pre_adapter_summary_metrics = step0001_picard_multiple_metrics.pre_adapter_summary_metrics
        File picard_pre_adapter_detail_metrics = step0001_picard_multiple_metrics.pre_adapter_detail_metrics

        Array[File] picard_wgs_metrics = step0002_picard_wgs_metrics.wgs_metrics

        File samtools_idxstats = step1001_samtools_idxstats.idxstats
        File samtools_flagstat = step1002_samtools_flagstat.flagstat
    }

}


task step0001_picard_multiple_metrics {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        File cram
        File cram_index
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
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        gatk CollectMultipleMetrics \
            --REFERENCE_SEQUENCE ~{reference_fasta} \
            --INPUT ~{cram} \
            --OUTPUT ~{output_prefix} \
            --PROGRAM null \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM QualityScoreDistribution \
            --PROGRAM MeanQualityByCycle \
            --PROGRAM CollectBaseDistributionByCycle \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectSequencingArtifactMetrics
    >>>

    output {
        File alignment_summary_metrics = "${output_prefix}.alignment_summary_metrics"
        File insert_size_metrics = "${output_prefix}.insert_size_metrics"
        File insert_size_histogram_pdf = "${output_prefix}.insert_size_histogram.pdf"
        File quality_distribution_metrics = "${output_prefix}.quality_distribution_metrics"
        File quality_distribution_pdf = "${output_prefix}.quality_distribution.pdf"
        File quality_by_cycle_metrics = "${output_prefix}.quality_by_cycle_metrics"
        File quality_by_cycle_pdf = "${output_prefix}.quality_by_cycle.pdf"
        File base_distribution_by_cycle_metrics = "${output_prefix}.base_distribution_by_cycle_metrics"
        File base_distribution_by_cycle_pdf = "${output_prefix}.base_distribution_by_cycle.pdf"
        File gc_bias_pdf = "${output_prefix}.gc_bias.pdf"
        File gc_bias_summary_metrics = "${output_prefix}.gc_bias.summary_metrics"
        File gc_bias_detail_metrics = "${output_prefix}.gc_bias.detail_metrics"
        File bait_bias_summary_metrics = "${output_prefix}.bait_bias_summary_metrics"
        File bait_bias_detail_metrics = "${output_prefix}.bait_bias_detail_metrics"
        File pre_adapter_summary_metrics = "${output_prefix}.pre_adapter_summary_metrics"
        File pre_adapter_detail_metrics = "${output_prefix}.pre_adapter_detail_metrics"
    }

}


task step0002_picard_wgs_metrics {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        String region_name
        File interval_list
        File cram
        File cram_index
        String wgs_metrics__name

        String docker_image
        Int threads
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

        base_name=$(echo ~{basename(cram)} | awk -F'.' '{ print $1 }')
        ln -s ~{cram} ${base_name}.~{region_name}.cram
        ln -s ~{cram_index} ${base_name}.~{region_name}.cram.crai

        gatk CollectWgsMetrics \
            --REFERENCE_SEQUENCE ~{reference_fasta} \
            --INPUT ${base_name}.~{region_name}.cram \
            --INTERVALS ~{interval_list} \
            --OUTPUT ~{wgs_metrics__name} \
    >>>

    output {
        File wgs_metrics = "${wgs_metrics__name}"
    }

}


task step1001_samtools_idxstats {

    input {
        File cram
        File cram_index
        String idxstats__name

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

        samtools idxstats \
            --verbosity 3 \
            --threads ~{threads} \
            ~{cram} \
        > ~{idxstats__name}
    >>>

    output {
        File idxstats = "${idxstats__name}"
    }

}


task step1002_samtools_flagstat {

    input {
        File cram
        File cram_index
        String flagstat__name

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

        samtools flagstat \
            --verbosity 3 \
            --threads ~{threads} \
            ~{cram} \
        > ~{flagstat__name}
    >>>

    output {
        File flagstat = "${flagstat__name}"
    }

}
