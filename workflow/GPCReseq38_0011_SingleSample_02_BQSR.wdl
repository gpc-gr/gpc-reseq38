#
#
#

version 1.0

import "./modules/md5sum.wdl"


workflow GPCReseq38_0011_SingleSample_02_BQSR {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes = [
            "${reference_fasta}.fai",
            sub(reference_fasta, ".fa(sta)?$", ".dict")
        ]

        String bqsr_model_region_type = "autosome"
        File bqsr_model_region_list = "${reference_fasta}.regions.${bqsr_model_region_type}.tsv"

        Array[VCF] known_site_vcfs

        String sample_id
        File sample_cram
        File sample_cram_index = "${sample_cram}.crai"

        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String base_recalibrator_docker_image = gatk_docker_image
        Int base_recalibrator_threads = 1
        Float base_recalibrator_memory_gb = 8
        String gather_bqsr_reports_docker_image = gatk_docker_image
        Int gather_bqsr_reports_threads = 1
        Float gather_bqsr_reports_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # preparation
    # --------------------------------------------------------------------------------

    Array[Object] regions = read_objects(bqsr_model_region_list)

    scatter (entry in known_site_vcfs) {
        File known_site_vcf = entry.vcf
        File known_site_vcf_index = entry.vcf_index
    }

    # --------------------------------------------------------------------------------
    # BQSR
    # --------------------------------------------------------------------------------

    scatter (region in regions) {

        call step1001_base_recalibrator { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            known_site_vcfs = known_site_vcf,
            known_site_vcf_indexes = known_site_vcf_index,
            region = "${region.contig}:${region.start}-${region.end}",
            cram = sample_cram,
            cram_index = sample_cram_index,
            bqsr_table__name = "${sample_id}.BQSR.table.${region.name}",
            docker_image = base_recalibrator_docker_image,
            threads = base_recalibrator_threads,
            memory_gb = base_recalibrator_memory_gb
        }

    }

    call step1002_gather_bqsr_reports { input:
        sources = step1001_base_recalibrator.bqsr_table,
        bqsr_table__name = "${sample_id}.BQSR.table.${bqsr_model_region_type}",
        docker_image = gather_bqsr_reports_docker_image,
        threads = gather_bqsr_reports_threads,
        memory_gb = gather_bqsr_reports_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call md5sum.md5sum as step9999_md5sum { input:
        sources = [step1002_gather_bqsr_reports.bqsr_table],
        md5sum_txt__name = "${sample_id}.GPCReseq_0011_SingleSample_02_BQSR.md5sum.txt"
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File bqsr_table = step1002_gather_bqsr_reports.bqsr_table
        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct VCF {
    File vcf
    File vcf_index
}


task step1001_base_recalibrator {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        Array[File] known_site_vcfs
        Array[File] known_site_vcf_indexes
        String region
        File cram
        File cram_index
        String bqsr_table__name

        String docker_image
        Int threads
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "${memory_gb}G"
        maxRetries: 3
    }

    command <<<
        set -euxo pipefail
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        gatk BaseRecalibrator \
            --reference ~{reference_fasta} \
            --input ~{cram} \
            --intervals ~{region} \
            ~{sep=" " prefix("--known-sites ", known_site_vcfs)} \
            --output ~{bqsr_table__name}
    >>>

    output {
        File bqsr_table = "${bqsr_table__name}"
    }

}


task step1002_gather_bqsr_reports {

    input {
        Array[File] sources
        String bqsr_table__name

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

        gatk GatherBQSRReports \
            ~{sep=" " prefix("--input ", sources)} \
            --output ~{bqsr_table__name}
    >>>

    output {
        File bqsr_table = "${bqsr_table__name}"
    }

}
