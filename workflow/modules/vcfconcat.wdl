#
#
#

version 1.0


task vcfconcat {

    input {
        Array[File] source_vcfs
        Array[File] source_vcf_indexes
        String result_vcf_gz__name
        Boolean apply_sort
        Boolean source_vcfs_have_same_header = false
        Boolean create_index = true

        String docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        Int threads = 4
        Float memory_gb = 1
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        apply_sort=~{true="yes" false="no" apply_sort}
        if [ "${apply_sort}" = "yes" ]; then
            bcftools concat \
                --no-version \
                --threads ~{threads} \
                --file-list ~{write_lines(source_vcfs)} \
                ~{true="--naive" false="" source_vcfs_have_same_header} \
                --output-type u \
            | bcftools sort \
                --temp-dir . \
                --output-type z \
                --output ~{result_vcf_gz__name}

        else
            bcftools concat \
                --no-version \
                --threads ~{threads} \
                --file-list ~{write_lines(source_vcfs)} \
                ~{true="--naive" false="" source_vcfs_have_same_header} \
                --output-type z \
                --output ~{result_vcf_gz__name}
        fi

        create_index=~{true="yes" false="no" create_index}
        if [ "${create_index}" = "yes" ]; then
            bcftools index \
                --threads ~{threads} \
                --tbi \
                ~{result_vcf_gz__name}
        else
            touch ~{result_vcf_gz__name}.tbi
        fi
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}
