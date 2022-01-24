#
#
#

version 1.0

import "./modules/cram2fastq.wdl"
import "./modules/cram2gvcf_auto.wdl"
import "./modules/fastq2cram_auto.wdl"
import "./modules/md5sum.wdl"


workflow GPCReseq38_0011_SingleSample_13_chrXY_PAR3 {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File baseline_reference_fasta
        Array[File] baseline_reference_fasta_general_indexes = [
            "${baseline_reference_fasta}.fai",
            sub(baseline_reference_fasta, ".fa(sta)?$", ".dict")
        ]

        File chrXY_PAR3_reference_fasta
        Array[File] chrXY_PAR3_reference_fasta_general_indexes = [
            "${chrXY_PAR3_reference_fasta}.fai",
            sub(chrXY_PAR3_reference_fasta, ".fa(sta)?$", ".dict")
        ]
        Array[File] chrXY_PAR3_reference_fasta_bwa2_indexes = [
            "${chrXY_PAR3_reference_fasta}.0123",
            "${chrXY_PAR3_reference_fasta}.amb",
            "${chrXY_PAR3_reference_fasta}.ann",
            "${chrXY_PAR3_reference_fasta}.bwt.2bit.64",
            "${chrXY_PAR3_reference_fasta}.pac",
            "${chrXY_PAR3_reference_fasta}.alt"
        ]

        String sample_id
        File sample_baseline_cram
        File sample_baseline_cram_index = "${sample_baseline_cram}.crai"
        File? sample_baseline_bqsr_table
        Array[String] ignored_read_groups = []

        String? accelaration

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
        Float bam2fastq_memory_gb = 16

        String align_cpu_docker_image = bwa_docker_image
        Int align_cpu_bwa_mem_threads = 8
        Int align_cpu_samtools_sort_threads = 4
        Float align_cpu_memory_gb = 32
        String rmdup_cpu_docker_image = gatk_docker_image
        Int rmdup_cpu_threads = 1
        Float rmdup_cpu_memory_gb = 8
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

    # --------------------------------------------------------------------------------
    # extract reads
    # --------------------------------------------------------------------------------

    call cram2fastq.cram2fastq as step0001_extract_reads { input:
        reference_fasta = baseline_reference_fasta,
        reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
        cram = sample_baseline_cram,
        cram_index = sample_baseline_cram_index,
        ignored_read_groups = ignored_read_groups,
        target_region_lists = [
            "${baseline_reference_fasta}.regions.chrXY_PAR2.chrX_nonPAR.tsv",
            "${baseline_reference_fasta}.regions.chrXY_PAR2.chrX_PAR.tsv",
            "${baseline_reference_fasta}.regions.chrXY_PAR2.chrY_nonPAR.tsv"
        ],
        python_docker_image = python_docker_image,
        extract_reads_docker_image = extract_reads_docker_image,
        extract_reads_threads = extract_reads_threads,
        extract_reads_memory_gb = extract_reads_memory_gb,
        bam2fastq_docker_image = bam2fastq_docker_image,
        bam2fastq_threads = bam2fastq_threads,
        bam2fastq_memory_gb = bam2fastq_memory_gb
    }

    # --------------------------------------------------------------------------------
    # alignment
    # --------------------------------------------------------------------------------

    call fastq2cram_auto.fastq2cram_auto as step0002_fastq2cram { input:
        reference_fasta = chrXY_PAR3_reference_fasta,
        reference_fasta_general_indexes = chrXY_PAR3_reference_fasta_general_indexes,
        reference_fasta_bwa2_indexes = chrXY_PAR3_reference_fasta_bwa2_indexes,
        sample_id = sample_id,
        read_pairs = step0001_extract_reads.read_pairs,
        cram__name = "${sample_id}.chrXY_PAR3.cram",
        accelaration = accelaration,
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

    # --------------------------------------------------------------------------------
    # variant call
    # --------------------------------------------------------------------------------

    String sample_id_with_suffix = sample_id + (if (defined(sample_baseline_bqsr_table)) then ".BQSR" else ".noBQSR")
    Array[VariantCallConfig] configs = [
        {
            "region_list": "${baseline_reference_fasta}.regions.chrXY_PAR3.chrX_nonPAR.tsv",
            "ploidy": 1,
            "gvcf_gz__name": "${sample_id_with_suffix}.chrXY_PAR3.chrX_nonPAR.male.g.vcf.gz"
        },
        {
            "region_list": "${baseline_reference_fasta}.regions.chrXY_PAR3.chrX_nonPAR.tsv",
            "ploidy": 2,
            "gvcf_gz__name": "${sample_id_with_suffix}.chrXY_PAR3.chrX_nonPAR.female.g.vcf.gz"
        },
        {
            "region_list": "${baseline_reference_fasta}.regions.chrXY_PAR3.chrX_PAR.tsv",
            "ploidy": 2,
            "gvcf_gz__name": "${sample_id_with_suffix}.chrXY_PAR3.chrX_PAR.g.vcf.gz"
        },
        {
            "region_list": "${baseline_reference_fasta}.regions.chrXY_PAR3.chrY_nonPAR.tsv",
            "ploidy": 1,
            "gvcf_gz__name": "${sample_id_with_suffix}.chrXY_PAR3.chrY_nonPAR.male.g.vcf.gz"
        }
    ]

    scatter (config in configs) {

        call cram2gvcf_auto.cram2gvcf_auto as step0003_cram2gvcf { input:
            reference_fasta = chrXY_PAR3_reference_fasta,
            reference_fasta_general_indexes = chrXY_PAR3_reference_fasta_general_indexes,
            regions = read_objects(config.region_list),
            sample_id = sample_id,
            sample_ploidy = config.ploidy,
            cram = step0002_fastq2cram.cram,
            cram_index = step0002_fastq2cram.cram_index,
            bqsr_table = sample_baseline_bqsr_table,
            gvcf_gz__name = config.gvcf_gz__name,
            accelaration = accelaration,
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

    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            step0003_cram2gvcf.gvcf_gz,
            step0003_cram2gvcf.gvcf_gz_index
        ]),
        md5sum_txt__name = "${sample_id_with_suffix}.GPCReseq_0011_SingleSample_13_chrXY_PAR3.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        Array[File] reads = step0001_extract_reads.reads

        File cram = step0002_fastq2cram.cram
        File cram_index = step0002_fastq2cram.cram_index

        Array[File] gvcf_gz = step0003_cram2gvcf.gvcf_gz
        Array[File] gvcf_gz_index = step0003_cram2gvcf.gvcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct VariantCallConfig {
    File region_list
    Int ploidy
    String gvcf_gz__name
}
