#
#
#

version 1.0

import "./modules/md5sum.wdl"
import "./modules/region2chunks.wdl"
import "./modules/vcfconcat.wdl"


workflow GPCReseq38_0013_Joint_09_ExtractSamplesFromVCF {

    input {
        String analysis_id
        String output_key

        Array[VCF] source_vcfs
        File sample_id_list
        Boolean force_samples

        Int chunk_size = 3 * 1000 * 1000    # 3Mb

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String python_docker_image = "python:3.8.6-slim-buster"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1
        String extract_samples_docker_image = bcftools_docker_image
        Int extract_samples_threads = 8
        Float extract_samples_memory_gb = 4
        String concat_vcfs_docker_image = bcftools_docker_image
        Int concat_vcfs_threads = 4
        Float concat_vcfs_memory_gb = 4
    }

    scatter (source in source_vcfs) {

        call region2chunks.region2chunks as step0000_split_region { input:
            contig = source.contig,
            start = source.start,
            end = source.end,
            chunk_size = chunk_size,
            chunks_tsv__name = "${basename(source.vcf)}.chunks.tsv",
            docker_image = split_region_docker_image,
            memory_gb = split_region_memory_gb
        }

        scatter (chunk in step0000_split_region.chunks) {

            call step0001_extract_samples { input:
                source_vcf = source.vcf,
                source_vcf_index = source.vcf_index,
                sample_id_list = sample_id_list,
                force_samples = force_samples,
                region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                result_vcf_gz__name = sub(basename(source.vcf), ".vcf(.gz)?", ".${output_key}.${chunk.name}.vcf.gz"),
                docker_image = extract_samples_docker_image,
                threads = extract_samples_threads,
                memory_gb = extract_samples_memory_gb
            }

        }

        call vcfconcat.vcfconcat as step0002_concat_vcfs { input:
            source_vcfs = step0001_extract_samples.result_vcf_gz,
            source_vcf_indexes = step0001_extract_samples.result_vcf_gz_index,
            result_vcf_gz__name = sub(basename(source.vcf), ".vcf(.gz)?", ".${output_key}.vcf.gz"),
            apply_sort = false,
            source_vcfs_have_same_header = true,
            create_index = true,
            docker_image = concat_vcfs_docker_image,
            threads = concat_vcfs_threads,
            memory_gb = concat_vcfs_memory_gb
        }

    }

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            step0002_concat_vcfs.result_vcf_gz,
            step0002_concat_vcfs.result_vcf_gz_index
        ]),
        md5sum_txt__name = "${analysis_id}.GPCReseq38_0013_Joint_09_ExtractSamplesFromVCF.md5sum.txt",
        docker_image = python_docker_image
    }

    output {
        Array[File] result_vcf_gz = step0002_concat_vcfs.result_vcf_gz
        Array[File] result_vcf_gz_index = step0002_concat_vcfs.result_vcf_gz_index

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct VCF {
    File vcf
    File vcf_index
    String contig
    Int start
    Int end
}


task step0001_extract_samples {

    input {
        File source_vcf
        File source_vcf_index
        File sample_id_list
        Boolean force_samples
        String region
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

        bcftools view \
            --no-version \
            --threads ~{threads} \
            --samples-file ~{sample_id_list} \
            ~{true="--force-samples" false="" force_samples} \
            --regions ~{region} \
            --targets ~{region} \
            --trim-alt-alleles \
            --output-type z \
            --output-file ~{result_vcf_gz__name}

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
