#
#
#

version 1.0


task report {

    input {
        Array[File] sources
        String output_prefix

        Boolean interactive = true
        Boolean generate_plots_zip = false

        String docker_image = "quay.io/biocontainers/multiqc:1.9--py_1"
        Int threads = 1
        Float memory_gb = 16
    }

    runtime {
        docker: docker_image
        cpu: threads
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        multiqc --file-list ~{write_lines(sources)} --fullnames ~{true="--interactive" false="" interactive} ~{true="--export" false="" generate_plots_zip}
        mv multiqc_report.html ~{output_prefix}.html

        python -c 'import shutil; shutil.make_archive("multiqc_data", "zip", "multiqc_data")'
        mv multiqc_data.zip ~{output_prefix}.zip

        generate_plots_zip=~{true="yes" false="no" generate_plots_zip}
        if [ "${generate_plots_zip}" = "yes" ]; then
            python -c 'import shutil; shutil.make_archive("multiqc_plots", "zip", "multiqc_plots")'
            mv multiqc_plots.zip ~{output_prefix}.plots.zip
        fi

        rm -rf multiqc_data
    >>>

    output {
        File html = "${output_prefix}.html"
        File zip = "${output_prefix}.zip"
        File? plots_zip = "${output_prefix}.plots.zip"
    }

}
