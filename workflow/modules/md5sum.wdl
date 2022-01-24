#
#
#

version 1.0


task md5sum {

    input {
        Array[File] sources
        String md5sum_txt__name

        String docker_image = "python:3.8.6-slim-buster"
        Float memory_gb = 1
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        md5sum $(cat ~{write_lines(sources)}) \
            | awk '{ l = split($2, path, "/"); print $1 "  " path[l] }' \
            | sort -k2,2 \
            > ~{md5sum_txt__name}
    >>>

    output {
        File md5sum_txt = "${md5sum_txt__name}"
    }

}
