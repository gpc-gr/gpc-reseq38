#
#
#

version 1.0


workflow GPCReseq38_0015_Phasing_01_Phasing {

    input {
        String analysis_id

        Array[VCF] source_vcfs
        File genetic_map_txt_gz

        Int window = 100000
        Int margin = 20000
        Float scaffold_maf_ge = 0.01
        Float scaffold_hwe_ge = 0.0001
        Float scaffold_missing_le = 0.05

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String python_docker_image = "python:3.8.6-slim-buster"
        String shapeit_docker_image = "quay.io/biocontainers/shapeit4:4.2.2--hc94963d_0"

        String prepare_chunks_docker_image = python_docker_image
        Int prepare_chunks_threads = 1
        Float prepare_chunks_memory_gb = 4
        String extract_sites_for_scaffold_docker_image = bcftools_docker_image
        Int extract_sites_for_scaffold_threads = 1
        Float extract_sites_for_scaffold_memory_gb = 4
        String phase_scaffold_docker_image = shapeit_docker_image
        Int phase_scaffold_threads = 16
        Float phase_scaffold_memory_gb = 32         # 64 -> 32
        String phase_with_scaffold_docker_image = shapeit_docker_image
        Int phase_with_scaffold_threads = 24
        Float phase_with_scaffold_memory_gb = 64    # 92 -> 64
        String ligate_docker_image = bcftools_docker_image
        Int ligate_threads = 4
        Float ligate_memory_gb = 16
    }

    scatter (vcf in source_vcfs) {

        call step0001_prepare_chunks { input:
            vcf = select_first([vcf.sites_vcf, vcf.vcf]),
            chunks_tsv__name = "${basename(vcf.vcf)}.chunks.tsv",
            window = window,
            margin = margin,
            docker_image = prepare_chunks_docker_image,
            threads = prepare_chunks_threads,
            memory_gb = prepare_chunks_memory_gb
        }

        scatter (chunk in step0001_prepare_chunks.chunks) {

            call step1001_extract_sites_for_scaffold { input:
                source_vcf = vcf.vcf,
                source_vcf_index = vcf.vcf_index,
                result_vcf_gz__name = sub(basename(vcf.vcf), ".vcf(.gz)?$", ".${chunk.name}.scaffold.unphased.vcf.gz"),
                region = "${chunk.contig}:${chunk.start_padded}-${chunk.end_padded}",
                maf_ge = scaffold_maf_ge,
                hwe_ge = scaffold_hwe_ge,
                missing_le = scaffold_missing_le,
                docker_image = extract_sites_for_scaffold_docker_image,
                threads = extract_sites_for_scaffold_threads,
                memory_gb = extract_sites_for_scaffold_memory_gb
            }

            call step1002_phase_scaffold { input:
                source_vcf = step1001_extract_sites_for_scaffold.result_vcf_gz,
                source_vcf_index = step1001_extract_sites_for_scaffold.result_vcf_gz_index,
                result_vcf_gz__name = sub(basename(vcf.vcf), ".vcf(.gz)?$", ".${chunk.name}.scaffold.phased.vcf.gz"),
                region = "${chunk.contig}:${chunk.start_padded}-${chunk.end_padded}",
                docker_image = phase_scaffold_docker_image,
                threads = phase_scaffold_threads,
                memory_gb = phase_scaffold_memory_gb
            }

            call step1003_phase_with_scaffold { input:
                source_vcf = vcf.vcf,
                source_vcf_index = vcf.vcf_index,
                scaffold_vcf = step1002_phase_scaffold.result_vcf_gz,
                scaffold_vcf_index = step1002_phase_scaffold.result_vcf_gz_index,
                result_vcf_gz__name = sub(basename(vcf.vcf), ".vcf(.gz)?$", ".${chunk.name}.phased.vcf.gz"),
                region = "${chunk.contig}:${chunk.start_padded}-${chunk.end_padded}",
                docker_image = phase_with_scaffold_docker_image,
                threads = phase_with_scaffold_threads,
                memory_gb = phase_with_scaffold_memory_gb
            }

        }

        call step2001_ligate { input:
            source_vcfs = step1003_phase_with_scaffold.result_vcf_gz,
            source_vcf_indexes = step1003_phase_with_scaffold.result_vcf_gz_index,
            result_vcf_gz__name = sub(basename(vcf.vcf), ".vcf(.gz)?$", ".phased.vcf.gz"),
            docker_image = ligate_docker_image,
            threads = ligate_threads,
            memory_gb = ligate_memory_gb
        }

    }

    output {
        Array[File] result_vcf_gz = step2001_ligate.result_vcf_gz
        Array[File] result_vcf_gz_index = step2001_ligate.result_vcf_gz_index
    }

}


struct VCF {
    File vcf
    File vcf_index
    File? sites_vcf
    File? sites_vcf_index
}


task step0001_prepare_chunks {

    input {
        File vcf
        String chunks_tsv__name
        Int window
        Int margin

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

        cat <<EOS >script.py
        if True:
            import sys

            #
            contig = None
            positions = []

            for line in sys.stdin:
                if line.startswith('#'):
                    continue

                cols = line.strip().split('\t')
                if contig:
                    assert cols[0] == contig
                else:
                    contig = cols[0]

                positions.append(cols[1])

            #
            window = ~{window}
            margin = ~{margin}
            index = 0

            print('\t'.join(['name', 'contig', 'start', 'end', 'start_padded', 'end_padded']))

            for chunk_start_index in range(0, len(positions), window):
                chunk_end_index = min(chunk_start_index + window, len(positions) - 1)
                padded_chunk_start_index = max(0, chunk_start_index - margin)
                padded_chunk_end_index = min(chunk_end_index + margin, len(positions) - 1)

                print('\t'.join([
                    f'chunk{index:06d}',
                    contig,
                    positions[chunk_start_index],
                    positions[chunk_end_index],
                    positions[padded_chunk_start_index],
                    positions[padded_chunk_end_index]
                ]))
                index += 1
        EOS

        gzip -dc ~{vcf} | python3 script.py > ~{chunks_tsv__name}
        rm script.py
    >>>

    output {
        File chunk_tsv = "${chunks_tsv__name}"
        Array[Object] chunks = read_objects(chunks_tsv__name)
    }

}


task step1001_extract_sites_for_scaffold {

    input {
        File source_vcf
        File source_vcf_index
        String result_vcf_gz__name

        String region
        Float maf_ge
        Float hwe_ge
        Float missing_le

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
            --regions ~{region} \
            --targets ~{region} \
            --apply-filters PASS \
            --types snps \
            --min-alleles 2 \
            --max-alleles 2 \
            --min-af ~{maf_ge}:minor \
            --output-type u \
            ~{source_vcf} \
        | bcftools +fill-tags \
            --no-version \
            --output-type u \
            -- \
            --tags HWE,F_MISSING \
        | bcftools view \
            --no-version \
            --threads ~{threads} \
            --include 'INFO/HWE >= ~{hwe_ge} && INFO/F_MISSING <= ~{missing_le}' \
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


task step1002_phase_scaffold {

    input {
        File source_vcf
        File source_vcf_index
        File genetic_map_txt_gz
        String result_vcf_gz__name

        String region

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

        shapeit4 \
            --thread ~{threads} \
            --seed 1 \
            --input ~{source_vcf} \
            --map ~{genetic_map_txt_gz} \
            --region ~{region} \
            --output ~{result_vcf_gz__name}

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


task step1003_phase_with_scaffold {

    input {
        File source_vcf
        File source_vcf_index
        File scaffold_vcf
        File scaffold_vcf_index
        File genetic_map_txt_gz
        String result_vcf_gz__name

        String region

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

        shapeit4 \
            --thread ~{threads} \
            --seed 1 \
            --sequencing \
            --input ~{source_vcf} \
            --map ~{genetic_map_txt_gz} \
            --region ~{region} \
            --scaffold ~{scaffold_vcf} \
            --output ~{result_vcf_gz__name}

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


task step2001_ligate {

    input {
        Array[File] source_vcfs
        Array[File] source_vcf_indexes
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

        bcftools concat \
            --no-version \
            --threads ~{threads} \
            --file-list ~{write_lines(source_vcfs)} \
            --ligate \
            --output-type z \
            --output ~{result_vcf_gz__name}

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
