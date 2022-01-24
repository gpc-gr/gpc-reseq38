#
#
#

version 1.0

import "./md5sum.wdl"
import "./vcfconcat.wdl"


workflow mitochondria_merge {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        File baseline_reference_fasta
        Array[File] baseline_reference_fasta_general_indexes = [
            "${baseline_reference_fasta}.fai",
            sub(baseline_reference_fasta, ".fa(sta)?$", ".dict")
        ]

        File mitochondria_shifted_reference_fasta
        Array[File] mitochondria_shifted_reference_fasta_general_indexes = [
            "${mitochondria_shifted_reference_fasta}.fai",
            sub(mitochondria_shifted_reference_fasta, ".fa(sta)?$", ".dict")
        ]

        Int mitochondria_size = read_objects("${baseline_reference_fasta}.regions.mitochondria.tsv")[0].end
        Int mitochondria_shift_size = read_int("${baseline_reference_fasta}.configs.mitochondria_shift_size.txt")

        File baseline_vcf
        File baseline_vcf_index
        File shifted_vcf
        File shifted_vcf_index

        String output_prefix

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String extract_snv_docker_image = bcftools_docker_image
        Int extract_snv_threads = 1
        Float extract_snv_memory_gb = 2

        String shift_back_docker_image = python_docker_image
        Int shift_back_threads = 2
        Float shift_back_memory_gb = 4

        String sort_vcf_docker_image = bcftools_docker_image
        Int sort_vcf_threads = 2
        Float sort_vcf_memory_gb = 4

        String compare_vcf_docker_image = bcftools_docker_image
        Int compare_vcf_threads = 4
        Float compare_vcf_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # baseline
    # --------------------------------------------------------------------------------

    call step0001_baseline_extract_snv { input:
        reference_fasta = baseline_reference_fasta,
        reference_fasta_general_indexes = baseline_reference_fasta_general_indexes,
        source_vcf = baseline_vcf,
        source_vcf_index = baseline_vcf_index,
        result_vcf_gz__name = "${output_prefix}.baseline.snv.vcf.gz",
        docker_image = extract_snv_docker_image,
        threads = extract_snv_threads,
        memory_gb = extract_snv_memory_gb
    }

    # --------------------------------------------------------------------------------
    # shifted
    # --------------------------------------------------------------------------------

    call step0001_baseline_extract_snv as step1001_shifted_extract_snv { input:
        reference_fasta = mitochondria_shifted_reference_fasta,
        reference_fasta_general_indexes = mitochondria_shifted_reference_fasta_general_indexes,
        source_vcf = shifted_vcf,
        source_vcf_index = shifted_vcf_index,
        result_vcf_gz__name = "${output_prefix}.shifted.snv.vcf.gz",
        docker_image = extract_snv_docker_image,
        threads = extract_snv_threads,
        memory_gb = extract_snv_memory_gb
    }

    call step1002_shifted_shift_back { input:
        mitochondria_size = mitochondria_size,
        mitochondria_shift_size = mitochondria_shift_size,
        source_vcf_gz = step1001_shifted_extract_snv.result_vcf_gz,
        result_vcf_gz__name = "${output_prefix}.shifted_back.snv.vcf.gz",
        docker_image = shift_back_docker_image,
        threads = shift_back_threads,
        memory_gb = shift_back_memory_gb
    }

    call step1003_shifted_sort_vcf { input:
        source_vcf_gz = step1002_shifted_shift_back.result_vcf_gz,
        result_vcf_gz__name = basename(step1002_shifted_shift_back.result_vcf_gz),
        docker_image = sort_vcf_docker_image,
        threads = sort_vcf_threads,
        memory_gb = sort_vcf_memory_gb
    }

    # --------------------------------------------------------------------------------
    # merge
    # --------------------------------------------------------------------------------

    call step2001_compare { input:
        baseline_vcf_gz = step0001_baseline_extract_snv.result_vcf_gz,
        baseline_vcf_gz_index = step0001_baseline_extract_snv.result_vcf_gz_index,
        shifted_back_vcf_gz = step1003_shifted_sort_vcf.result_vcf_gz,
        shifted_back_vcf_gz_index = step1003_shifted_sort_vcf.result_vcf_gz_index,
        output_prefix = output_prefix,
        docker_image = compare_vcf_docker_image,
        threads = compare_vcf_threads,
        memory_gb = compare_vcf_memory_gb
    }

    call vcfconcat.vcfconcat as step2002_merge { input:
        source_vcfs = [
            step2001_compare.common_baseline_vcf_gz,
            step2001_compare.unique_baseline_vcf_gz,
            step2001_compare.unique_shifted_back_vcf_gz
        ],
        source_vcf_indexes = [
            step2001_compare.common_baseline_vcf_gz_index,
            step2001_compare.unique_baseline_vcf_gz_index,
            step2001_compare.unique_shifted_back_vcf_gz_index
        ],
        result_vcf_gz__name = "${output_prefix}.merged",
        apply_sort = true,
        source_vcfs_have_same_header = false,
        create_index = true
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File merged_vcf_gz = step2002_merge.result_vcf_gz
        File merged_vcf_gz_index = step2002_merge.result_vcf_gz_index
        File common_baseline_vcf_gz = step2001_compare.common_baseline_vcf_gz
        File common_baseline_vcf_gz_index = step2001_compare.common_baseline_vcf_gz_index
        File common_shifted_back_vcf_gz = step2001_compare.common_shifted_back_vcf_gz
        File common_shifted_back_vcf_gz_index = step2001_compare.common_shifted_back_vcf_gz_index
        File unique_baseline_vcf_gz = step2001_compare.unique_baseline_vcf_gz
        File unique_baseline_vcf_gz_index = step2001_compare.unique_baseline_vcf_gz_index
        File unique_shifted_back_vcf_gz = step2001_compare.unique_shifted_back_vcf_gz
        File unique_shifted_back_vcf_gz_index = step2001_compare.unique_shifted_back_vcf_gz_index
    }

}


task step0001_baseline_extract_snv {

    input {
        File reference_fasta
        Array[File] reference_fasta_general_indexes

        File source_vcf
        File source_vcf_index
        String result_vcf_gz__name

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

        bcftools norm \
            --no-version \
            --threads 1 \
            --fasta-ref ~{reference_fasta} \
            --multiallelics -any \
            --output-type u \
            ~{source_vcf} \
        | bcftools view \
            --no-version \
            --threads 1 \
            --types snps \
            --output-type u \
        | bcftools norm \
            --no-version \
            --threads 1 \
            --fasta-ref ~{reference_fasta} \
            --multiallelics -any \
            --output-type z \
            --output ~{result_vcf_gz__name}

        bcftools index \
            --threads 1 \
            --tbi \
            ~{result_vcf_gz__name}
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}


task step1002_shifted_shift_back {

    input {
        Int mitochondria_size
        Int mitochondria_shift_size

        File source_vcf_gz
        String result_vcf_gz__name

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

        #
        cat <<CODE >fix-position.py
        if True:
            import sys

            MITOCHONDRIA_SIZE = ~{mitochondria_size}
            mitochondria_shift_size = ~{mitochondria_shift_size}

            for line in sys.stdin:
                if not line.startswith('#'):
                    cols = list(line.split('\t'))

                    position = int(cols[1]) + mitochondria_shift_size
                    if position > MITOCHONDRIA_SIZE:
                        position -= MITOCHONDRIA_SIZE

                    cols[1] = str(position)
                    line = '\t'.join(cols)

                sys.stdout.write(line)
        CODE

        #
        gzip -dc ~{source_vcf_gz} | python3 fix-position.py | gzip -c > ~{result_vcf_gz__name}
        rm fix-position.py
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
    }

}


task step1003_shifted_sort_vcf {

    input {
        File source_vcf_gz
        String result_vcf_gz__name

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

        bcftools sort \
            --max-mem ~{floor(memory_gb * 1000 * 8/10)}M \
            --temp-dir . \
            --output-type z \
            --output ~{result_vcf_gz__name} \
            ~{source_vcf_gz}

        bcftools index \
            --threads ~{threads} \
            --tbi \
            ~{result_vcf_gz__name}
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}


task step2001_compare {

    input {
        File baseline_vcf_gz
        File baseline_vcf_gz_index
        File shifted_back_vcf_gz
        File shifted_back_vcf_gz_index
        String output_prefix

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
        # compare VCFs
        mkdir -p ./outputs
        bcftools isec \
            --no-version \
            --threads ~{threads} \
            --prefix ./outputs \
            ~{baseline_vcf_gz} \
            ~{shifted_back_vcf_gz}

        # compress results of bcftools isec
        function compress_and_index () {
            bcftools view --no-version --threads ~{threads} --output-type z --output-file $2 $1
            bcftools index --threads ~{threads} --tbi $2
        }

        compress_and_index ./outputs/0000.vcf ~{output_prefix}.unique.baseline.vcf.gz
        compress_and_index ./outputs/0001.vcf ~{output_prefix}.unique.shifted_back.vcf.gz
        compress_and_index ./outputs/0002.vcf ~{output_prefix}.common.baseline.vcf.gz
        compress_and_index ./outputs/0003.vcf ~{output_prefix}.common.shifted_back.vcf.gz
    >>>

    output {
        File common_baseline_vcf_gz = "${output_prefix}.common.baseline.vcf.gz"
        File common_baseline_vcf_gz_index = "${output_prefix}.common.baseline.vcf.gz.tbi"
        File common_shifted_back_vcf_gz = "${output_prefix}.common.shifted_back.vcf.gz"
        File common_shifted_back_vcf_gz_index = "${output_prefix}.common.shifted_back.vcf.gz.tbi"
        File unique_baseline_vcf_gz = "${output_prefix}.unique.baseline.vcf.gz"
        File unique_baseline_vcf_gz_index = "${output_prefix}.unique.baseline.vcf.gz.tbi"
        File unique_shifted_back_vcf_gz = "${output_prefix}.unique.shifted_back.vcf.gz"
        File unique_shifted_back_vcf_gz_index = "${output_prefix}.unique.shifted_back.vcf.gz.tbi"
    }

}
