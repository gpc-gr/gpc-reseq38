#
#
#

version 1.0


task report {

    input {
        File read

        String docker_image = "quay.io/biocontainers/fastqc:0.11.9--0"
        Int threads = 4
        Float memory_gb = 4
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        fastqc \
            --threads ~{threads} \
            --nogroup \
            --outdir . \
            --dir . \
            ~{read}
    >>>

    output {
        File zip = sub(basename(read), ".f(ast)?q.gz$", "_fastqc.zip")
    }

}
