#
#
#

version 1.0


workflow cram2gvcf_cpu {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        Array[Object] regions

        String sample_id
        Int sample_ploidy
        File cram
        File cram_index
        File? bqsr_table
        Float? snv_heterozygosity
        Float? indel_heterozygosity

        String gvcf_gz__name

        String apply_bqsr_cpu_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int apply_bqsr_cpu_threads = 1
        Float apply_bqsr_cpu_memory_gb = 8
        String haplotype_caller_cpu_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int haplotype_caller_cpu_threads = 4
        Float haplotype_caller_cpu_memory_gb = 4
        String concat_gvcf_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        Int concat_gvcf_threads = 4
        Float concat_gvcf_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # variant call
    # --------------------------------------------------------------------------------

    scatter (region in regions) {

        if (defined(bqsr_table)) {

            call step0001_apply_bqsr { input:
                reference_fasta = reference_fasta,
                reference_fasta_general_indexes = reference_fasta_general_indexes,
                raw_cram = cram,
                raw_cram_index = cram_index,
                bqsr_table = select_first([bqsr_table]),
                region = "${region.contig}:${region.start}-${region.end}",
                bqsr_cram__name = "${sample_id}.BQSR.${region.name}.cram",
                docker_image = apply_bqsr_cpu_docker_image,
                threads = apply_bqsr_cpu_threads,
                memory_gb = apply_bqsr_cpu_memory_gb
            }

        }

        File cram_for_variant_call = select_first([step0001_apply_bqsr.bqsr_cram, cram])
        File cram_index_for_variant_call = select_first([step0001_apply_bqsr.bqsr_cram_index, cram_index])
        String chunk_gvcf__name = "${sample_id}" + (if (defined(bqsr_table)) then ".BQSR." else ".noBQSR." ) + "${region.name}.g.vcf"

        call step1001_haplotype_caller { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            region = "${region.contig}:${region.start}-${region.end}",
            cram = cram_for_variant_call,
            cram_index = cram_index_for_variant_call,
            ploidy = sample_ploidy,
            snv_heterozygosity = snv_heterozygosity,
            indel_heterozygosity = indel_heterozygosity,
            gvcf__name = chunk_gvcf__name,
            docker_image = haplotype_caller_cpu_docker_image,
            threads = haplotype_caller_cpu_threads,
            memory_gb = haplotype_caller_cpu_memory_gb
        }

    }

    call step1002_concat_gvcfs { input:
        sources = step1001_haplotype_caller.gvcf,
        source_indexes = step1001_haplotype_caller.gvcf_index,
        gvcf_gz__name = gvcf_gz__name,
        docker_image = concat_gvcf_docker_image,
        threads = concat_gvcf_threads,
        memory_gb = concat_gvcf_memory_gb
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File gvcf_gz = step1002_concat_gvcfs.gvcf_gz
        File gvcf_gz_index = step1002_concat_gvcfs.gvcf_gz_index
    }

}


task step0001_apply_bqsr {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        File raw_cram
        File raw_cram_index
        File bqsr_table
        String region
        String bqsr_cram__name

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

        gatk ApplyBQSR \
            --reference ~{reference_fasta} \
            --input ~{raw_cram} \
            --intervals ~{region} \
            --use-original-qualities \
            --bqsr-recal-file ~{bqsr_table} \
            --static-quantized-quals 10 \
            --static-quantized-quals 20 \
            --static-quantized-quals 30 \
            --output ~{bqsr_cram__name}
    >>>

    output {
        File bqsr_cram = "${bqsr_cram__name}"
        File bqsr_cram_index = "${bqsr_cram__name}.bai"
    }

}


task step1001_haplotype_caller {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        String region
        File cram
        File cram_index
        Int ploidy
        Float? snv_heterozygosity
        Float? indel_heterozygosity
        String gvcf__name

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

        gatk HaplotypeCaller \
            --reference ~{reference_fasta} \
            --input ~{cram} \
            --sample-ploidy ~{ploidy} \
            ~{"--heterozygosity " + snv_heterozygosity} \
            ~{"--indel-heterozygosity " + indel_heterozygosity} \
            --intervals ~{region} \
            --output ~{gvcf__name} \
            --emit-ref-confidence GVCF \
            --native-pair-hmm-threads ~{threads}
    >>>

    output {
        File gvcf = "${gvcf__name}"
        File gvcf_index = "${gvcf__name}.idx"
    }

}


task step1002_concat_gvcfs {

    input {
        Array[File] sources
        Array[File] source_indexes

        String gvcf_gz__name

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

        bcftools concat \
            --no-version \
            --threads ~{threads} \
            --file-list ~{write_lines(sources)} \
            --output-type z \
            --output ~{gvcf_gz__name}

        bcftools index \
            --threads ~{threads} \
            --tbi \
            ~{gvcf_gz__name}
    >>>

    output {
        File gvcf_gz = "${gvcf_gz__name}"
        File gvcf_gz_index = "${gvcf_gz__name}.tbi"
    }

}
