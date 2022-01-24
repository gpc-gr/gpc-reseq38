#
#
#

version 1.0

import "./modules/multiqc.wdl"


workflow SeqQC_01_FastQC {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        String flowcell
        Array[ReadPair] read_pairs

        String? accelaration

        String falco_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-falco:falco_0.2.4"
        Int falco_threads = 1
        Float falco_memory_gb = 1

        String fastqc_docker_image = "quay.io/biocontainers/fastqc:0.11.9--0"
        Int fastqc_threads = 1
        Float fastqc_memory_gb = 2

        String multiqc_docker_image = "quay.io/biocontainers/multiqc:1.9--py_1"
        Int multiqc_threads = 2
        Float multiqc_memory_gb = 48
    }

    # --------------------------------------------------------------------------------
    # FastQC
    # --------------------------------------------------------------------------------

    scatter (read_pair in read_pairs) {

        scatter (read in [read_pair.read1, read_pair.read2]) {

            if (defined(accelaration) && accelaration == "FALCO") {

                call step0001_falco { input:
                    read = read,
                    docker_image = falco_docker_image,
                    threads = falco_threads,
                    memory_gb = falco_memory_gb
                }

            }

            if (!defined(accelaration)) {

                call step0001_fastqc { input:
                    read = read,
                    docker_image = fastqc_docker_image,
                    threads = fastqc_threads,
                    memory_gb = fastqc_memory_gb
                }

            }

        }

    }

    Array[File] fastqc_htmls = select_all(flatten(flatten([step0001_falco.html, step0001_fastqc.html])))
    Array[File] fastqc_zips = select_all(flatten(flatten([step0001_falco.zip, step0001_fastqc.zip])))

    # --------------------------------------------------------------------------------
    # MultiQC
    # --------------------------------------------------------------------------------

    call multiqc.report as step1001_multiqc { input:
        sources = fastqc_zips,
        output_prefix = "${flowcell}.01_FastQC.multiqc",
        generate_plots_zip = true,
        threads = multiqc_threads,
        memory_gb = multiqc_memory_gb
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File multiqc_html = step1001_multiqc.html
        File multiqc_zip = step1001_multiqc.zip
        File multiqc_plots_zip = step1001_multiqc.plots_zip
    }

}


struct ReadPair {
    File read1
    File read2
}


task step0001_falco {

    input {
        File read

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

        falco \
            --nogroup \
            --outdir . \
            --dir . \
            ~{read}

        falco-make-inline-html \
            fastqc_report.html \
            ~{sub(basename(read), ".f(ast)?q.gz$", "_fastqc.html")}

        falco-make-zip \
            --fastq-name ~{basename(read)} \
            --fastqc-data-txt fastqc_data.txt \
            --fastqc-report-html fastqc_report.html \
            --summary-txt summary.txt \
            --output-directory .
    >>>

    output {
        File html = sub(basename(read), ".f(ast)?q.gz$", "_fastqc.html")
        File zip = sub(basename(read), ".f(ast)?q.gz$", "_fastqc.zip")
    }

}


task step0001_fastqc {

    input {
        File read

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

        fastqc \
            --threads ~{threads} \
            --nogroup \
            --outdir . \
            --dir . \
            ~{read}
    >>>

    output {
        File html = sub(basename(read), ".f(ast)?q.gz$", "_fastqc.html")
        File zip = sub(basename(read), ".f(ast)?q.gz$", "_fastqc.zip")
    }

}
