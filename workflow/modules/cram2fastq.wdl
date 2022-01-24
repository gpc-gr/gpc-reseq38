#
#
#

version 1.0


workflow cram2fastq {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        File cram
        File cram_index
        Array[String] ignored_read_groups

        Array[File] target_region_lists

        String python_docker_image = "python:3.8.6-slim-buster"
        String extract_reads_docker_image = "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
        Int extract_reads_threads = 6
        Float extract_reads_memory_gb = 4
        String bam2fastq_docker_image = "broadinstitute/gatk:4.1.0.0"
        Int bam2fastq_threads = 1
        Float bam2fastq_memory_gb = 8
    }

    call step0001_collect_target_contigs { input:
        cram = cram,
        target_region_lists = target_region_lists,
        docker_image = python_docker_image
    }

    call step1001_extract_reads_from_cram { input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        cram = cram,
        cram_index = cram_index,
        target_contigs = step0001_collect_target_contigs.target_contigs,
        bam__name = "${basename(cram)}.extracted.bam",
        bam_header_txt__name = "${basename(cram)}.extracted.bam.header.txt",
        docker_image = extract_reads_docker_image,
        threads = extract_reads_threads,
        memory_gb = extract_reads_memory_gb
    }

    call step1002_convert_bam_to_fastq { input:
        bam = step1001_extract_reads_from_cram.bam,
        bam_header_txt = step1001_extract_reads_from_cram.bam_header_txt,
        ignored_read_groups = ignored_read_groups,
        docker_image = bam2fastq_docker_image,
        threads = bam2fastq_threads,
        memory_gb = bam2fastq_memory_gb
    }

    output {
        Array[File] reads = step1002_convert_bam_to_fastq.reads
        Array[Pair[Pair[String, String], Pair[File, File]]] read_pairs = step1002_convert_bam_to_fastq.read_pairs
    }

}


task step0001_collect_target_contigs {

    input {
        File cram
        Array[File] target_region_lists

        String docker_image
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "1G"
    }

    command <<<
        set -euxo pipefail

        cat ~{write_lines(target_region_lists)} \
            | xargs -n1 awk 'NR > 1 { print $2 }' \
            | sort \
            | uniq \
            > contigs.txt
    >>>

    output {
        Array[String] target_contigs = read_lines("contigs.txt")
    }

}


task step1001_extract_reads_from_cram {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes
        File cram
        File cram_index
        Array[String] target_contigs
        String bam__name
        String bam_header_txt__name

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

        # extract reads
        cram0=~{cram}
        bam1=${TMPDIR}/unmapped1.bam
        bam2=${TMPDIR}/unmapped2.bam
        bam3=${TMPDIR}/unmapped3.bam
        bam4=${TMPDIR}/mapped.bam

        samtools view --threads ~{threads} -T ~{reference_fasta} -b -u -f  4 -F 264 -o ${bam1} ${cram0}
        samtools view --threads ~{threads} -T ~{reference_fasta} -b -u -f  8 -F 260 -o ${bam2} ${cram0}
        samtools view --threads ~{threads} -T ~{reference_fasta} -b -u -f 12 -F 256 -o ${bam3} ${cram0}
        samtools view --threads ~{threads} -T ~{reference_fasta} -b -u       -F 268 -o ${bam4} ${cram0} ~{sep=" " target_contigs}

        # merge BAMs
        samtools merge \
            --threads ~{ceil(threads / 2)} \
            -c \
            -p \
            - \
            ${bam1} ${bam2} ${bam3} ${bam4} \
        | samtools sort \
            --threads ~{ceil(threads / 2)} \
            -T ${TMPDIR}/$(basename ${cram0}).sorting \
            -n \
            -o ~{bam__name} \
            -

        #
        samtools view -H ${cram0} > ~{bam_header_txt__name}
    >>>

    output {
        File bam = "${bam__name}"
        File bam_header_txt = "${bam_header_txt__name}"
    }

}


task step1002_convert_bam_to_fastq {

    input {
        File bam
        File bam_header_txt
        Array[String] ignored_read_groups

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

        # convert BAM to FASTQs
        gatk SamToFastq \
            --INPUT ~{bam} \
            --OUTPUT_PER_RG \
            --OUTPUT_DIR $(pwd) \
            --INCLUDE_NON_PF_READS \
            --RG_TAG ID \
            --VALIDATION_STRINGENCY SILENT

        # collect RGs
        cat ~{bam_header_txt} | grep ^@RG | sed 's/\t/\\t/g' | sort > outputs.RG.raw.txt

        # collect FASTQs
        find . -name "*.fastq" | xargs -n1 -I%p sh -c 'gzip %p'

        cat <<EOS | python3 -
        if True:
            import collections
            import glob
            import gzip
            import os

            fastqs = collections.defaultdict(list)
            for path in sorted(glob.glob('./*.fastq.gz')):
                #
                ok = False
                with gzip.open(path, 'rt') as fin:
                    for line in fin:
                        if line.strip():
                            ok = True
                            break

                if not ok:
                    continue

                #
                id = os.path.basename(path).replace('_1.fastq.gz', '').replace('_2.fastq.gz', '')
                fastqs[id].append(path)

            fout_id = open('outputs.ID.txt', 'w')
            fout_rg = open('outputs.RG.txt', 'w')
            fout_r1 = open('outputs.R1.txt', 'w')
            fout_r2 = open('outputs.R2.txt', 'w')
            fout_files = open('outputs.FILES.txt', 'w')

            for line in open('outputs.RG.raw.txt'):
                line = line.strip()
                if not line:
                    continue

                id = [e[3:] for e in line.split('\\\\t') if e.startswith('ID:')][0]
                if (id in fastqs) and (len(fastqs[id]) == 2):
                    print(id, file=fout_id)
                    print(line, file=fout_rg)

                    read1, read2 = list(sorted(fastqs[id]))
                    print(read1, file=fout_r1)
                    print(read2, file=fout_r2)
                    print(read1, file=fout_files)
                    print(read2, file=fout_files)
        EOS
    >>>

    output {
        Array[File] reads = read_lines("outputs.FILES.txt")
        Array[Pair[Pair[String, String], Pair[File, File]]] read_pairs = zip(
            zip(
                read_lines("outputs.ID.txt"),
                read_lines("outputs.RG.txt")
            ),
            zip(
                read_lines("outputs.R1.txt"),
                read_lines("outputs.R2.txt")
            )
        )
    }

}
