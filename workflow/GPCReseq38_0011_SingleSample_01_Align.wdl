#
#
#

version 1.0

import "./modules/cramqc.wdl"
import "./modules/fastq2cram_auto.wdl"
import "./modules/md5sum.wdl"
import "./modules/multiqc.wdl"
import "./modules/types.wdl"


workflow GPCReseq38_0011_SingleSample_01_Align {

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

        Array[String] wgs_metrics_regions = [
            "autosome",
            "chrXY_PAR2.chrX_nonPAR",
            "chrXY_PAR2.chrY_nonPAR",
            "mitochondria"
        ]

        String sample_id
        Array[ReadPair] sample_read_pairs

        String? accelaration

        String bwa_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-bwamem2:bwamem2_2.2.1-samtools_1.11"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String multiqc_docker_image = "quay.io/biocontainers/multiqc:1.9--py_1"
        String python_docker_image = "python:3.8.6-slim-buster"
        String samtools_docker_image = "quay.io/biocontainers/samtools:1.11--h6270b1f_0"

        Int multiqc_threads = 1
        Float multiqc_memory_gb = 16

        String align_cpu_docker_image = bwa_docker_image
        Int align_cpu_total_threads = 10
        Int align_cpu_bwa_mem_threads = 8
        Int align_cpu_samtools_sort_threads = 4
        Float align_cpu_memory_gb = 32
        String rmdup_cpu_docker_image = gatk_docker_image
        Int rmdup_cpu_threads = 1
        Float rmdup_cpu_memory_gb = 16
        String bam2cram_docker_image = samtools_docker_image
        Int bam2cram_threads = 8
        Float bam2cram_memory_gb = 1

        String picard_multiple_metrics_docker_image = gatk_docker_image
        Int picard_multiple_metrics_threads = 1
        Float picard_multiple_metrics_memory_gb = 8
        String picard_wgs_metrics_docker_image = gatk_docker_image
        Int picard_wgs_metrics_threads = 1
        Float picard_wgs_metrics_memory_gb = 4
        String samtools_idxstats_docker_image = samtools_docker_image
        Int samtools_idxstats_threads = 8
        Float samtools_idxstats_memory_gb = 8
        String samtools_flagstat_docker_image = samtools_docker_image
        Int samtools_flagstat_threads = 8
        Float samtools_flagstat_memory_gb = 8
    }

    # --------------------------------------------------------------------------------
    # alignment
    # --------------------------------------------------------------------------------

    call types.convert_read_pairs_to_internal_representation as step1001_convert_read_pairs_to_internal_representation { input:
        sample_id = sample_id,
        read_pairs = sample_read_pairs,
        docker_image = python_docker_image
    }

    call fastq2cram_auto.fastq2cram_auto as step1002_fastq2cram { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        reference_fasta_bwa2_indexes = reference_fasta_bwa2_indexes,
        sample_id = sample_id,
        read_pairs = step1001_convert_read_pairs_to_internal_representation.read_pairs_internal,
        cram__name = "${sample_id}.base.cram",
        accelaration = accelaration,
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

    # --------------------------------------------------------------------------------
    # alignment QC
    # --------------------------------------------------------------------------------

    call cramqc.cramqc as step1003_cramqc { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        wgs_metrics_regions = wgs_metrics_regions,
        cram = step1002_fastq2cram.cram,
        cram_index = step1002_fastq2cram.cram_index,
        picard_multiple_metrics_docker_image = picard_multiple_metrics_docker_image,
        picard_multiple_metrics_threads = picard_multiple_metrics_threads,
        picard_multiple_metrics_memory_gb = picard_multiple_metrics_memory_gb,
        picard_wgs_metrics_docker_image = picard_wgs_metrics_docker_image,
        picard_wgs_metrics_threads = picard_wgs_metrics_threads,
        picard_wgs_metrics_memory_gb = picard_wgs_metrics_memory_gb,
        samtools_idxstats_docker_image = samtools_idxstats_docker_image,
        samtools_idxstats_threads = samtools_idxstats_threads,
        samtools_idxstats_memory_gb = samtools_idxstats_memory_gb,
        samtools_flagstat_docker_image = samtools_flagstat_docker_image,
        samtools_flagstat_threads = samtools_flagstat_threads,
        samtools_flagstat_memory_gb = samtools_flagstat_memory_gb
    }

    Array[File] base_cram_metrics = flatten([
        [
            step1002_fastq2cram.rmdup_metrics,
            step1003_cramqc.picard_alignment_summary_metrics,
            step1003_cramqc.picard_insert_size_metrics,
            step1003_cramqc.picard_quality_distribution_metrics,
            step1003_cramqc.picard_quality_by_cycle_metrics,
            step1003_cramqc.picard_base_distribution_by_cycle_metrics,
            step1003_cramqc.picard_gc_bias_summary_metrics,
            step1003_cramqc.picard_gc_bias_detail_metrics,
            step1003_cramqc.picard_bait_bias_summary_metrics,
            step1003_cramqc.picard_bait_bias_detail_metrics,
            step1003_cramqc.picard_pre_adapter_summary_metrics,
            step1003_cramqc.picard_pre_adapter_detail_metrics,
            step1003_cramqc.samtools_idxstats,
            step1003_cramqc.samtools_flagstat
        ],
        step1003_cramqc.picard_wgs_metrics
    ])
    Array[File] base_cram_metrics_figures = [
        step1003_cramqc.picard_insert_size_histogram_pdf,
        step1003_cramqc.picard_quality_distribution_pdf,
        step1003_cramqc.picard_quality_by_cycle_pdf,
        step1003_cramqc.picard_base_distribution_by_cycle_pdf,
        step1003_cramqc.picard_gc_bias_pdf
    ]

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call multiqc.report as step9001_multiqc { input:
        sources = base_cram_metrics,
        output_prefix = "${sample_id}.base.multiqc",
        docker_image = multiqc_docker_image,
        threads = multiqc_threads,
        memory_gb = multiqc_memory_gb
    }

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            base_cram_metrics,
            [
                step1002_fastq2cram.cram,
                step1002_fastq2cram.cram_index,
                step9001_multiqc.html,
                step9001_multiqc.zip
            ]
        ]),
        md5sum_txt__name = "${sample_id}.GPCReseq_0011_SingleSample_01_Align.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File cram = step1002_fastq2cram.cram
        File cram_index = step1002_fastq2cram.cram_index
        Array[File] cram_metrics = base_cram_metrics
        Array[File] cram_metrics_figures = base_cram_metrics_figures

        File multiqc_html = step9001_multiqc.html
        File multiqc_zip = step9001_multiqc.zip

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}
