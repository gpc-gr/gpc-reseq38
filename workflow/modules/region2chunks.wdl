#
#
#

version 1.0


task region2chunks {

    input {
        String contig
        Int start
        Int end
        Int chunk_size
        String chunks_tsv__name

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

        python3 - <<EOS >~{chunks_tsv__name}
        if True:
            contig = '~{contig}'
            start = ~{start}
            end = ~{end}
            size = ~{chunk_size}
            index = 1

            print('\t'.join(['name', 'contig', 'start', 'end']))
            while start <= end:
                print('\t'.join(map(str, [
                    'chunk{:08d}'.format(index),
                    contig,
                    start,
                    min(start + size - 1, end)
                ])))

                start += size
                index += 1
        EOS
    >>>

    output {
        File chunks_tsv = "${chunks_tsv__name}"
        Array[Object] chunks = read_objects(chunks_tsv__name)
    }

}

