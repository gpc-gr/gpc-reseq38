#
#
#

version 1.0

import "./modules/cram2fastq.wdl"
import "./modules/cram2gvcf_auto.wdl"
import "./modules/fastq2cram_auto.wdl"
import "./modules/md5sum.wdl"


workflow GPCReseq38_0011_SingleSample_14_Mitochondria {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File baseline_reference_fasta
        Array[File] baseline_reference_fasta_bwa2_indexes = [
            "${baseline_reference_fasta}.0123",
            "${baseline_reference_fasta}.amb",
            "${baseline_reference_fasta}.ann",
            "${baseline_reference_fasta}.bwt.2bit.64",
            "${baseline_reference_fasta}.pac",
            "${baseline_reference_fasta}.alt"
        ]
        Array[File] baseline_reference_fasta_general_indexes = [
            "${baseline_reference_fasta}.fai",
            sub(baseline_reference_fasta, ".fa(sta)?$", ".dict")
        ]

        File mitochondria_shifted_reference_fasta
        Array[File] mitochondria_shifted_reference_fasta_general_indexes = [
            "${mitochondria_shifted_reference_fasta}.fai",
            sub(mitochondria_shifted_reference_fasta, ".fa(sta)?$", ".dict")
        ]
        Array[File] mitochondria_shifted_reference_fasta_bwa2_indexes = [
            "${mitochondria_shifted_reference_fasta}.0123",
            "${mitochondria_shifted_reference_fasta}.amb",
            "${mitochondria_shifted_reference_fasta}.ann",
            "${mitochondria_shifted_reference_fasta}.bwt.2bit.64",
            "${mitochondria_shifted_reference_fasta}.pac",
            "${mitochondria_shifted_reference_fasta}.alt"
        ]

        String sample_id
        File sample_baseline_cram
        File sample_baseline_cram_index
        File? sample_baseline_bqsr_table
        Array[String] ignored_read_groups = []
        Float snv_heterozygosity = 0.0025
        Float indel_heterozygosity = 0.0003125

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String bwa_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-bwamem2:bwamem2_2.2.1-samtools_1.11"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"
        String samtools_docker_image = "quay.io/biocontainers/samtools:1.11--h6270b1f_0"

        String extract_reads_docker_image = samtools_docker_image
        Int extract_reads_threads = 6
        Float extract_reads_memory_gb = 4
        String bam2fastq_docker_image = gatk_docker_image
        Int bam2fastq_threads = 1
        Float bam2fastq_memory_gb = 8

        String align_cpu_docker_image = bwa_docker_image
        Int align_cpu_bwa_mem_threads = 8
        Int align_cpu_samtools_sort_threads = 4
        Float align_cpu_memory_gb = 32
        String rmdup_cpu_docker_image = gatk_docker_image
        Int rmdup_cpu_threads = 1
        Float rmdup_cpu_memory_gb = 16
        String bam2cram_docker_image = samtools_docker_image
        Int bam2cram_threads = 8
        Float bam2cram_memory_gb = 8

        String apply_bqsr_cpu_docker_image = gatk_docker_image
        Int apply_bqsr_cpu_threads = 1
        Float apply_bqsr_cpu_memory_gb = 8
        String haplotype_caller_cpu_docker_image = gatk_docker_image
        Int haplotype_caller_cpu_threads = 4
        Float haplotype_caller_cpu_memory_gb = 4
        String concat_gvcf_docker_image = bcftools_docker_image
        Int concat_gvcf_threads = 4
        Float concat_gvcf_memory_gb = 4
    }

    String sample_id_with_suffix = sample_id + (if (defined(sample_baseline_bqsr_table)) then ".BQSR" else ".noBQSR")

    # --------------------------------------------------------------------------------
    # extract reads
    # --------------------------------------------------------------------------------

    call cram2fastq.cram2fastq as step0001_extract_reads { input:
        reference_fasta = baseline_reference_fasta,
        reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
        cram = sample_baseline_cram,
        cram_index = sample_baseline_cram_index,
        ignored_read_groups = ignored_read_groups,
        target_region_lists = ["${baseline_reference_fasta}.regions.mitochondria.tsv"],
        python_docker_image = python_docker_image,
        extract_reads_docker_image = extract_reads_docker_image,
        extract_reads_threads = extract_reads_threads,
        extract_reads_memory_gb = extract_reads_memory_gb,
        bam2fastq_docker_image = bam2fastq_docker_image,
        bam2fastq_threads = bam2fastq_threads,
        bam2fastq_memory_gb = bam2fastq_memory_gb
    }

    # --------------------------------------------------------------------------------
    # alignment & variant call (baseline)
    # --------------------------------------------------------------------------------

    call fastq2cram_auto.fastq2cram_auto as step1001_fastq2cram_base { input:
        reference_fasta = baseline_reference_fasta,
        reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
        reference_fasta_bwa2_indexes = baseline_reference_fasta_bwa2_indexes,
        sample_id = sample_id,
        read_pairs = step0001_extract_reads.read_pairs,
        cram__name = "${sample_id}.mitochondria.base.cram",
        align_cpu_docker_image = align_cpu_docker_image,
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

    call cram2gvcf_auto.cram2gvcf_auto as step1002_cram2gvcf_base { input:
        reference_fasta = baseline_reference_fasta,
        reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
        regions = read_objects("${baseline_reference_fasta}.regions.mitochondria.tsv"),
        sample_id = sample_id,
        sample_ploidy = 1,
        cram = step1001_fastq2cram_base.cram,
        cram_index = step1001_fastq2cram_base.cram_index,
        bqsr_table = sample_baseline_bqsr_table,
        snv_heterozygosity = snv_heterozygosity,
        indel_heterozygosity = indel_heterozygosity,
        gvcf_gz__name = "${sample_id_with_suffix}.mitochondria.base.g.vcf.gz",
        apply_bqsr_cpu_docker_image = apply_bqsr_cpu_docker_image,
        apply_bqsr_cpu_threads = apply_bqsr_cpu_threads,
        apply_bqsr_cpu_memory_gb = apply_bqsr_cpu_memory_gb,
        haplotype_caller_cpu_docker_image = haplotype_caller_cpu_docker_image,
        haplotype_caller_cpu_threads = haplotype_caller_cpu_threads,
        haplotype_caller_cpu_memory_gb = haplotype_caller_cpu_memory_gb,
        concat_gvcf_docker_image = concat_gvcf_docker_image,
        concat_gvcf_threads = concat_gvcf_threads,
        concat_gvcf_memory_gb = concat_gvcf_memory_gb
    }

    # --------------------------------------------------------------------------------
    # alignment & variant call (shifted)
    # --------------------------------------------------------------------------------

    call fastq2cram_auto.fastq2cram_auto as step2001_fastq2cram_shifted { input:
        reference_fasta = mitochondria_shifted_reference_fasta,
        reference_fasta_general_indexes = mitochondria_shifted_reference_fasta_general_indexes,
        reference_fasta_bwa2_indexes = mitochondria_shifted_reference_fasta_bwa2_indexes,
        sample_id = sample_id,
        read_pairs = step0001_extract_reads.read_pairs,
        cram__name = "${sample_id}.mitochondria.shifted.cram",
        align_cpu_docker_image = align_cpu_docker_image,
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

    call cram2gvcf_auto.cram2gvcf_auto as step2002_cram2gvcf_shifted { input:
        reference_fasta = mitochondria_shifted_reference_fasta,
        reference_fasta_general_indexes = mitochondria_shifted_reference_fasta_general_indexes,
        regions = read_objects("${baseline_reference_fasta}.regions.mitochondria.tsv"),
        sample_id = sample_id,
        sample_ploidy = 1,
        cram = step2001_fastq2cram_shifted.cram,
        cram_index = step2001_fastq2cram_shifted.cram_index,
        bqsr_table = sample_baseline_bqsr_table,
        snv_heterozygosity = snv_heterozygosity,
        indel_heterozygosity = indel_heterozygosity,
        gvcf_gz__name = "${sample_id_with_suffix}.mitochondria.shifted.g.vcf.gz",
        apply_bqsr_cpu_docker_image = apply_bqsr_cpu_docker_image,
        apply_bqsr_cpu_threads = apply_bqsr_cpu_threads,
        apply_bqsr_cpu_memory_gb = apply_bqsr_cpu_memory_gb,
        haplotype_caller_cpu_docker_image = haplotype_caller_cpu_docker_image,
        haplotype_caller_cpu_threads = haplotype_caller_cpu_threads,
        haplotype_caller_cpu_memory_gb = haplotype_caller_cpu_memory_gb,
        concat_gvcf_docker_image = concat_gvcf_docker_image,
        concat_gvcf_threads = concat_gvcf_threads,
        concat_gvcf_memory_gb = concat_gvcf_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            step0001_extract_reads.reads,
            [
                step1001_fastq2cram_base.cram,
                step1001_fastq2cram_base.cram_index,
                step1002_cram2gvcf_base.gvcf_gz,
                step1002_cram2gvcf_base.gvcf_gz_index,
                step2001_fastq2cram_shifted.cram,
                step2001_fastq2cram_shifted.cram_index,
                step2002_cram2gvcf_shifted.gvcf_gz,
                step2002_cram2gvcf_shifted.gvcf_gz_index
            ]
        ]),
        md5sum_txt__name = "${sample_id}.GPCReseq_0011_SingleSample_14_Mitochondria.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        Array[File] reads = step0001_extract_reads.reads

        File base_cram = step1001_fastq2cram_base.cram
        File base_cram_index = step1001_fastq2cram_base.cram_index
        File base_gvcf_gz = step1002_cram2gvcf_base.gvcf_gz
        File base_gvcf_gz_index = step1002_cram2gvcf_base.gvcf_gz_index

        File shifted_cram = step2001_fastq2cram_shifted.cram
        File shifted_cram_index = step2001_fastq2cram_shifted.cram_index
        File shifted_gvcf_gz = step2002_cram2gvcf_shifted.gvcf_gz
        File shifted_gvcf_gz_index = step2002_cram2gvcf_shifted.gvcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}
