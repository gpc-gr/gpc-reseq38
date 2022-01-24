#
#
#

version 1.0


struct ReadPair {
    String id
    String library
    String platform_unit
    String platform
    String platform_model
    String sequencing_center_name
    File read1
    File read2
}


task convert_read_pairs_to_internal_representation {

    input {
        String sample_id
        Array[ReadPair] read_pairs

        String docker_image = "python:3.8.6-slim-buster"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "1G"
    }

    command <<<
        python3 <<CODE >fastqs.txt
        if True:
            import csv

            print('\t'.join(['id', 'rg', 'read1', 'read2']))

            with open('~{write_objects(read_pairs)}') as fin:
                for record in csv.DictReader(fin, delimiter='\t'):
                    print('\t'.join([
                        record['id'],
                        '@RG\\\\t' + '\\\\t'.join([
                            'ID:{id}'.format(**record),
                            'SM:~{sample_id}',
                            'PU:{platform_unit}'.format(**record),
                            'LB:{library}'.format(**record),
                            'PL:{platform}'.format(**record),
                            'PM:{platform_model}'.format(**record),
                            'CN:{sequencing_center_name}'.format(**record)
                        ]),
                        record['read1'],
                        record['read2']
                    ]))
        CODE

        cat fastqs.txt | awk 'NR > 1 { print $1 }' > fastqs.ID.txt
        cat fastqs.txt | awk 'NR > 1 { print $2 }' > fastqs.RG.txt
        cat fastqs.txt | awk 'NR > 1 { print $3 }' > fastqs.R1.txt
        cat fastqs.txt | awk 'NR > 1 { print $4 }' > fastqs.R2.txt
    >>>

    output {
        Array[Pair[Pair[String, String], Pair[File, File]]] read_pairs_internal = zip(
            zip(
                read_lines("fastqs.ID.txt"),
                read_lines("fastqs.RG.txt")
            ),
            zip(
                read_lines("fastqs.R1.txt"),
                read_lines("fastqs.R2.txt")
            )
        )
    }

}
