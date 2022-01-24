#
#
#

version 1.0

import "./md5sum.wdl"


workflow gvcf2vcf_direct {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes = [
            "${reference_fasta}.fai",
            sub(reference_fasta, ".fa(sta)?$", ".dict")
        ]

        File gvcf
        File gvcf_index
        String vcf_gz__name

        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String genotype_gvcfs_docker_image = gatk_docker_image
        Int genotype_gvcfs_threads = 2
        Float genotype_gvcfs_memory_gb = 4
    }

    call step0001_genotype_gvcfs { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        gvcf = gvcf,
        gvcf_index = gvcf_index,
        vcf_gz__name = vcf_gz__name,
        docker_image = genotype_gvcfs_docker_image,
        threads = genotype_gvcfs_threads,
        memory_gb = genotype_gvcfs_memory_gb
    }

    output {
        File vcf_gz = step0001_genotype_gvcfs.vcf_gz
        File vcf_gz_index = step0001_genotype_gvcfs.vcf_gz_index
    }

}


task step0001_genotype_gvcfs {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        File gvcf
        File gvcf_index

        String vcf_gz__name

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

        export JAVA_TOOL_OPTIONS="-XX:+UseSerialGC -XX:CICompilerCount=2 -Xmx~{floor(memory_gb * 1000 * 8/10)}m"
        gatk GenotypeGVCFs \
            --reference ~{reference_fasta} \
            --variant ~{gvcf} \
            --output ~{vcf_gz__name} \
            --annotation-group StandardAnnotation \
            --annotation-group AS_StandardAnnotation
    >>>

    output {
        File vcf_gz = "${vcf_gz__name}"
        File vcf_gz_index = "${vcf_gz__name}.tbi"
    }

}
