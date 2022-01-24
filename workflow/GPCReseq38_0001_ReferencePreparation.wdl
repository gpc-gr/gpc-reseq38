#
#
#

version 1.0

import "./modules/md5sum.wdl"


workflow GPCReseq38_0001_ReferencePreparation {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File reference_fasta
        File? reference_fasta_alt
        Boolean generate_baseline_indexes = true

        String chrY_contig_name
        Int chrY_XTR_start
        Int chrY_XTR_end

        String mitochondria_contig_name
        Int mitochondria_shift_size

        String bedtools_docker_image = "quay.io/biocontainers/bedtools:2.27.1--he513fc3_4"
        String bwa_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-bwamem2:bwamem2_2.2.1-samtools_1.11"
        String python_docker_image = "python:3.8.6-slim-buster"
        String samtools_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-samtools:perl_5.32.0-samtools_1.11"

        String bwa_index_docker_image = bwa_docker_image
        Int bwa_index_threads = 1
        Float bwa_index_memory_gb = 64
        String samtools_faidx_docker_image = samtools_docker_image
        Int samtools_faidx_threads = 1
        Float samtools_faidx_memory_gb = 8
        String samtools_dict_docker_image = samtools_docker_image
        Int samtools_dict_threads = 1
        Float samtools_dict_memory_gb = 8
        Float samtools_seq_cache_populate_memory_gb = 4
        String maskfasta_docker_image = bedtools_docker_image
        Int maskfasta_threads = 1
        Float maskfasta_memory_gb = 4
        String shift_sequence_docker_image = python_docker_image
        Int shift_sequence_threads = 1
        Float shift_sequence_memory_gb = 4
    }

    String reference_fasta_name_wo_extension = sub(basename(reference_fasta), ".fa(sta)?$", "")
    String reference_fasta_extension = sub(basename(reference_fasta), reference_fasta_name_wo_extension, "")

    # --------------------------------------------------------------------------------
    # baseline
    # --------------------------------------------------------------------------------

    if (generate_baseline_indexes) {

        call bwa_index as step0001_baseline_bwa_index { input:
            reference_fasta = reference_fasta,
            reference_fasta_alt = reference_fasta_alt,
            docker_image = bwa_index_docker_image,
            threads = bwa_index_threads,
            memory_gb = bwa_index_memory_gb
        }

        call samtools_faidx as step0002_baseline_samtools_faidx { input:
            reference_fasta = reference_fasta,
            docker_image = samtools_faidx_docker_image,
            threads = samtools_faidx_threads,
            memory_gb = samtools_faidx_memory_gb
        }

        call samtools_dict as step0003_baseline_samtools_dict { input:
            reference_fasta = reference_fasta,
            docker_image = samtools_dict_docker_image,
            threads = samtools_dict_threads,
            memory_gb = samtools_dict_memory_gb
        }

        call samtools_seq_cache_populate as step0004_baseline_samtools_seq_cache_populate { input:
            reference_fasta = reference_fasta,
            docker_image = samtools_dict_docker_image,
            memory_gb = samtools_seq_cache_populate_memory_gb
        }

    }

    # --------------------------------------------------------------------------------
    # chrXY_PAR3
    # --------------------------------------------------------------------------------

    call step1000_chrXY_PAR3_mask_XTR_on_chrY { input:
        reference_fasta = reference_fasta,
        chrY_contig_name = chrY_contig_name,
        chrY_XTR_start = chrY_XTR_start,
        chrY_XTR_end = chrY_XTR_end,
        masked_reference_fasta__name = "${reference_fasta_name_wo_extension}.chrXY_PAR3${reference_fasta_extension}",
        docker_image = maskfasta_docker_image,
        threads = maskfasta_threads,
        memory_gb = maskfasta_memory_gb
    }

    call bwa_index as step1001_chrXY_PAR3_bwa_index { input:
        reference_fasta = step1000_chrXY_PAR3_mask_XTR_on_chrY.masked_reference_fasta,
        reference_fasta_alt = reference_fasta_alt,
        docker_image = bwa_index_docker_image,
        threads = bwa_index_threads,
        memory_gb = bwa_index_memory_gb
    }

    call samtools_faidx as step1002_chrXY_PAR3_samtools_faidx { input:
        reference_fasta = step1000_chrXY_PAR3_mask_XTR_on_chrY.masked_reference_fasta,
        docker_image = samtools_faidx_docker_image,
        threads = samtools_faidx_threads,
        memory_gb = samtools_faidx_memory_gb
    }

    call samtools_dict as step1003_chrXY_PAR3_samtools_dict { input:
        reference_fasta = step1000_chrXY_PAR3_mask_XTR_on_chrY.masked_reference_fasta,
        docker_image = samtools_dict_docker_image,
        threads = samtools_dict_threads,
        memory_gb = samtools_dict_memory_gb
    }

    # --------------------------------------------------------------------------------
    # Mitochondria
    # --------------------------------------------------------------------------------

    call step2000_mitochondria_shift_mitochondria { input:
        reference_fasta = reference_fasta,
        mitochondria_contig_name = mitochondria_contig_name,
        mitochondria_shift_size = mitochondria_shift_size,
        mitochondria_shifted_reference_fasta__name = "${reference_fasta_name_wo_extension}.mitochondria.shifted${reference_fasta_extension}",
        docker_image = shift_sequence_docker_image,
        threads = shift_sequence_threads,
        memory_gb = shift_sequence_memory_gb
    }

    call bwa_index as step2001_mitochondria_bwa_index { input:
        reference_fasta = step2000_mitochondria_shift_mitochondria.mitochondria_shifted_reference_fasta,
        reference_fasta_alt = reference_fasta_alt,
        docker_image = bwa_index_docker_image,
        threads = bwa_index_threads,
        memory_gb = bwa_index_memory_gb
    }

    call samtools_faidx as step2002_mitochondria_samtools_faidx { input:
        reference_fasta = step2000_mitochondria_shift_mitochondria.mitochondria_shifted_reference_fasta,
        docker_image = samtools_faidx_docker_image,
        threads = samtools_faidx_threads,
        memory_gb = samtools_faidx_memory_gb
    }

    call samtools_dict as step2003_mitochondria_samtools_dict { input:
        reference_fasta = step2000_mitochondria_shift_mitochondria.mitochondria_shifted_reference_fasta,
        docker_image = samtools_dict_docker_image,
        threads = samtools_dict_threads,
        memory_gb = samtools_dict_memory_gb
    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    Array[File] baseline_files = if (generate_baseline_indexes)
        then select_all(flatten([
            [
                reference_fasta,
                step0002_baseline_samtools_faidx.reference_fasta_fai,
                step0003_baseline_samtools_dict.reference_dict,
                step0004_baseline_samtools_seq_cache_populate.reference_fasta_cache_tar_gz
            ],
            select_first([step0001_baseline_bwa_index.reference_fasta_bwa2_indexes])
        ]))
        else []

    Array[File] chrXY_PAR3_files = flatten([
        [
            step1000_chrXY_PAR3_mask_XTR_on_chrY.masked_reference_fasta,
            step1002_chrXY_PAR3_samtools_faidx.reference_fasta_fai,
            step1003_chrXY_PAR3_samtools_dict.reference_dict
        ],
        step1001_chrXY_PAR3_bwa_index.reference_fasta_bwa2_indexes
    ])

    Array[File] mitochondria_files = flatten([
        [
            step2000_mitochondria_shift_mitochondria.mitochondria_shifted_reference_fasta,
            step2002_mitochondria_samtools_faidx.reference_fasta_fai,
            step2003_mitochondria_samtools_dict.reference_dict
        ],
        step2001_mitochondria_bwa_index.reference_fasta_bwa2_indexes
    ])

    Array[File] all_files = select_all(flatten([
        baseline_files,
        chrXY_PAR3_files,
        mitochondria_files
    ]))

    call md5sum.md5sum as step9999_md5sum { input:
        sources = all_files,
        md5sum_txt__name = "${basename(reference_fasta)}.GPCReseq38_0001_ReferencePreparation.md5sum.txt"
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        Array[File] artifacts = all_files
        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


task bwa_index {

    input {
        File reference_fasta
        File? reference_fasta_alt

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

        bwa-mem2 index \
            -p ~{basename(reference_fasta)} \
            ~{reference_fasta}

        have_alt=~{true="yes" false="no" defined(reference_fasta_alt)}
        if [ "${have_alt}" = "yes" ]; then
            cp ~{reference_fasta_alt} ~{basename(reference_fasta)}.64.alt
        fi
    >>>

    output {
        Array[File?] reference_fasta_bwa2_indexes_with_nullable = [
            "${basename(reference_fasta)}.0123",
            "${basename(reference_fasta)}.amb",
            "${basename(reference_fasta)}.ann",
            "${basename(reference_fasta)}.bwt.2bit.64",
            "${basename(reference_fasta)}.pac",
            "${basename(reference_fasta)}.alt"
        ]
        Array[File] reference_fasta_bwa2_indexes = select_all(reference_fasta_bwa2_indexes_with_nullable)
    }

}


task samtools_faidx {

    input {
        File reference_fasta

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

        ln -s ~{reference_fasta} ~{basename(reference_fasta)}
        samtools faidx ~{basename(reference_fasta)}
    >>>

    output {
        File reference_fasta_fai = "${basename(reference_fasta)}.fai"
    }

}


task samtools_dict {

    input {
        File reference_fasta
        String reference_dict__name = sub(basename(reference_fasta), ".fa(sta)?$", ".dict")

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
        samtools dict ~{reference_fasta} > ~{reference_dict__name}
    >>>

    output {
        File reference_dict = "${reference_dict__name}"
    }

}


task samtools_seq_cache_populate {

    input {
        File reference_fasta
        String reference_fasta_cache_tar_gz__name = "${basename(reference_fasta)}.samtools.cache.tar.gz"

        String docker_image
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "${memory_gb}G"
    }

    command <<<
        set -euxo pipefail

        mkdir -p ./cache
        seq_cache_populate.pl -root ./cache ~{reference_fasta}

        tar zcvf ~{reference_fasta_cache_tar_gz__name} cache
    >>>

    output {
        File reference_fasta_cache_tar_gz = "${reference_fasta_cache_tar_gz__name}"
    }

}


task step1000_chrXY_PAR3_mask_XTR_on_chrY {

    input {
        File reference_fasta
        String chrY_contig_name
        Int chrY_XTR_start
        Int chrY_XTR_end
        String masked_reference_fasta__name

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

        bedtools maskfasta \
            -fi ~{reference_fasta} \
            -bed <(echo -e '~{chrY_contig_name}\t~{chrY_XTR_start-1}\t~{chrY_XTR_end-1}') \
            -fo ~{masked_reference_fasta__name}
    >>>

    output {
        File masked_reference_fasta = "${masked_reference_fasta__name}"
    }

}


task step2000_mitochondria_shift_mitochondria {

    input {
        File reference_fasta
        String mitochondria_contig_name
        Int mitochondria_shift_size
        String mitochondria_shifted_reference_fasta__name

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

        python3 <<CODE >~{mitochondria_shifted_reference_fasta__name}
        if True:
            # --------------------------------------------------------------------------------
            def parse_fasta(fin):
                header = None
                sequence = None

                for line in fin:
                    line = line.strip()
                    if not line:
                        continue

                    if line.startswith('>'):
                        if header:
                            yield header, sequence

                        header = line
                        sequence = []

                    else:
                        sequence.append(line)

                if header:
                    yield header, sequence

            # --------------------------------------------------------------------------------
            with open('~{reference_fasta}') as fin:
                for header, sequence in parse_fasta(fin):
                    if header[1:].split()[0] == '~{mitochondria_contig_name}':
                       original_sequence = ''.join(sequence)
                       shift_size = ~{mitochondria_shift_size}
                       shifted_sequence = original_sequence[shift_size:] + original_sequence[:shift_size]

                       width = len(sequence[0])
                       sequence = [shifted_sequence[s:s+width] for s in range(0, len(shifted_sequence), width)]

                    print(header)
                    for line in sequence:
                        print(line)

                    print()
        CODE
    >>>

    output {
        File mitochondria_shifted_reference_fasta = "${mitochondria_shifted_reference_fasta__name}"
    }

}
