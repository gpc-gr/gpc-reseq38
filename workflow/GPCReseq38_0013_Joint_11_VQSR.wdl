#
#
#

version 1.0

import "./modules/md5sum.wdl"
import "./modules/multiqc.wdl"
import "./modules/region2chunks.wdl"
import "./modules/vcfconcat.wdl"

#
# based on https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
# and https://github.com/gatk-workflows/broad-prod-wgs-germline-snps-indels/tree/1.2.0
#

workflow GPCReseq38_0013_Joint_11_VQSR {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        String analysis_id

        Array[ResourceVCF] snv_resource_vcfs
        Array[ResourceVCF] indel_resource_vcfs

        Array[Float] snv_tranches = [100.0, 99.95, 99.9, 99.8, 99.6, 99.5, 99.4, 99.3, 99.0, 98.0, 97.0, 90.0]
        Array[String] snv_annotations = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]
        Int snv_max_gaussians = 6
        Float snv_filter_level = 99.7
        Array[Float] indel_tranches = [100.0, 99.95, 99.9, 99.5, 99.0, 97.0, 96.0, 95.0, 94.0, 93.5, 93.0, 92.0, 91.0, 90.0]
        Array[String] indel_annotations = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]
        Int indel_max_gaussians = 4
        Float indel_filter_level = 99.7

        Array[InputVCF] source_vcfs
        Int chunk_size = 3 * 1000 * 1000    # 3Mb

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String gatk_docker_image = "broadinstitute/gatk:4.1.0.0"
        String multiqc_docker_image = "quay.io/biocontainers/multiqc:1.9--py_1"
        String python_docker_image = "python:3.8.6-slim-buster"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1
        String filter_by_excess_het_docker_image = bcftools_docker_image
        Int filter_by_excess_het_threads = 4
        Float filter_by_excess_het_memory_gb = 1
        String concat_vcfs_docker_image = bcftools_docker_image
        Int concat_vcfs_threads = 4
        Float concat_vcfs_memory_gb = 1
        String variant_recalibrator_docker_image = gatk_docker_image
        Float variant_recalibrator_memory_gb = 32
        String apply_vqsr_docker_image = gatk_docker_image
        Float apply_vqsr_memory_gb = 8
        String transfer_filtering_result_docker_image = bcftools_docker_image
        Int transfer_filtering_result_threads = 4
        Float transfer_filtering_result_memory_gb = 4
        String bcftools_stats_docker_image = bcftools_docker_image
        Int bcftools_stats_threads = 4
        Float bcftools_stats_memory_gb = 2
        String combine_bcftools_stats_docker_image = python_docker_image
        Int combine_bcftools_stats_threads = 1
        Float combine_bcftools_stats_memory_gb = 4
        Int multiqc_threads = 1
        Float multiqc_memory_gb = 16
    }

    # --------------------------------------------------------------------------------
    # preparation
    # --------------------------------------------------------------------------------

    scatter (entry in source_vcfs) {

        call region2chunks.region2chunks as step0000_split_region { input:
            contig = entry.contig,
            start = entry.start,
            end = entry.end,
            chunk_size = chunk_size,
            chunks_tsv__name = "${basename(entry.vcf)}.chunks.tsv",
            docker_image = split_region_docker_image,
            memory_gb = split_region_memory_gb
        }

        scatter (chunk in step0000_split_region.chunks) {

            call step0001_filter_by_excess_het { input:
                source_vcf = entry.vcf,
                source_vcf_index = entry.vcf_index,
                region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                result_vcf_gz__name = sub(basename(entry.vcf), ".vcf(.gz)?$", ".${chunk.name}.ExessHet.vcf.gz"),
                docker_image = filter_by_excess_het_docker_image,
                threads = filter_by_excess_het_threads,
                memory_gb = filter_by_excess_het_memory_gb
            }

        }

        call vcfconcat.vcfconcat as step0002_concat_vcfs { input:
            source_vcfs = step0001_filter_by_excess_het.result_vcf_gz,
            source_vcf_indexes = step0001_filter_by_excess_het.result_vcf_gz,
            result_vcf_gz__name = sub(basename(entry.vcf), ".vcf(.gz)?$", ".ExessHet.vcf.gz"),
            apply_sort = false,
            source_vcfs_have_same_header = true,
            create_index = true,
            docker_image = concat_vcfs_docker_image,
            threads = concat_vcfs_threads,
            memory_gb = concat_vcfs_memory_gb
        }

    }

    # --------------------------------------------------------------------------------
    # VQSR model construction
    # --------------------------------------------------------------------------------

    call step1000_variant_recalibrator as step1001_variant_recalibrator_indel { input:
        vcfs = step0002_concat_vcfs.result_vcf_gz,
        vcf_indexes = step0002_concat_vcfs.result_vcf_gz_index,
        tranches = indel_tranches,
        annotations = indel_annotations,
        mode = "INDEL",
        max_gaussians = indel_max_gaussians,
        resources = indel_resource_vcfs,
        output_prefix = "${analysis_id}.VQSR.INDEL",
        docker_image = variant_recalibrator_docker_image,
        memory_gb = variant_recalibrator_memory_gb
    }

    call step1000_variant_recalibrator as step1002_variant_recalibrator_snv { input:
        vcfs = step0002_concat_vcfs.result_vcf_gz,
        vcf_indexes = step0002_concat_vcfs.result_vcf_gz_index,
        tranches = snv_tranches,
        annotations = snv_annotations,
        mode = "SNP",
        max_gaussians = snv_max_gaussians,
        resources = snv_resource_vcfs,
        output_prefix = "${analysis_id}.VQSR.SNV",
        docker_image = variant_recalibrator_docker_image,
        memory_gb = variant_recalibrator_memory_gb
    }

    # --------------------------------------------------------------------------------
    # VQSR & stats
    # --------------------------------------------------------------------------------

    scatter (index in range(length(source_vcfs))) {

        InputVCF source_vcf = source_vcfs[index]
        File excess_het_vcf = step0002_concat_vcfs.result_vcf_gz[index]
        File excess_het_vcf_index = step0002_concat_vcfs.result_vcf_gz_index[index]

        call step2000_apply_vqsr as step2001_apply_vqsr_indel { input:
            source_vcf = excess_het_vcf,
            source_vcf_index = excess_het_vcf_index,
            recal_file = step1001_variant_recalibrator_indel.recal_file,
            recal_index_file = step1001_variant_recalibrator_indel.recal_index_file,
            tranches_file = step1001_variant_recalibrator_indel.tranches_file,
            filter_level = indel_filter_level,
            mode = "INDEL",
            result_vcf_gz__name = sub(basename(excess_het_vcf), ".vcf.gz$", ".VQSR.indel_only.vcf.gz"),
            docker_image = apply_vqsr_docker_image,
            memory_gb = apply_vqsr_memory_gb
        }

        call step2000_apply_vqsr as step2002_apply_vqsr_snv { input:
            source_vcf = step2001_apply_vqsr_indel.result_vcf_gz,
            source_vcf_index = step2001_apply_vqsr_indel.result_vcf_gz_index,
            recal_file = step1002_variant_recalibrator_snv.recal_file,
            recal_index_file = step1002_variant_recalibrator_snv.recal_index_file,
            tranches_file = step1002_variant_recalibrator_snv.tranches_file,
            filter_level = snv_filter_level,
            mode = "SNP",
            result_vcf_gz__name = sub(basename(excess_het_vcf), ".vcf.gz$", ".VQSR.vcf.gz"),
            docker_image = apply_vqsr_docker_image,
            memory_gb = apply_vqsr_memory_gb
        }

        scatter (chunk in step0000_split_region.chunks[index]) {

            call step2003_transfer_filtering_result { input:
                original_vcf = source_vcf.vcf,
                original_vcf_index = source_vcf.vcf_index,
                annotation_vcf = step2002_apply_vqsr_snv.result_vcf_gz,
                annotation_vcf_index = step2002_apply_vqsr_snv.result_vcf_gz_index,
                region = "${chunk.contig}:${chunk.start}-${chunk.end}",
                result_vcf_gz__name = sub(basename(source_vcf.vcf), ".vcf(.gz)?", "") + ".${chunk.name}.ExcessHet.VQSR.vcf.gz",
                docker_image = transfer_filtering_result_docker_image,
                threads = transfer_filtering_result_threads,
                memory_gb = transfer_filtering_result_memory_gb
            }

            call step2100_bcftools_stats as step2101_bcftools_stats_all { input:
                source_vcf = step2003_transfer_filtering_result.result_vcf_gz,
                source_vcf_index = step2003_transfer_filtering_result.result_vcf_gz_index,
                stats__name = "${basename(step2003_transfer_filtering_result.result_vcf_gz)}.bcftools.stats.ALL",
                docker_image = bcftools_stats_docker_image,
                threads = bcftools_stats_threads,
                memory_gb = bcftools_stats_memory_gb
            }

            call step2100_bcftools_stats as step2102_bcftools_stats_pass { input:
                source_vcf = step2003_transfer_filtering_result.result_vcf_gz,
                source_vcf_index = step2003_transfer_filtering_result.result_vcf_gz_index,
                filter = "PASS",
                stats__name = "${basename(step2003_transfer_filtering_result.result_vcf_gz)}.bcftools.stats.PASS",
                docker_image = bcftools_stats_docker_image,
                threads = bcftools_stats_threads,
                memory_gb = bcftools_stats_memory_gb
            }

        }

        call vcfconcat.vcfconcat as step2004_concat_vcfs { input:
            source_vcfs = step2003_transfer_filtering_result.result_vcf_gz,
            source_vcf_indexes = step2003_transfer_filtering_result.result_vcf_gz_index,
            result_vcf_gz__name = sub(basename(source_vcf.vcf), ".vcf(.gz)?", "") + ".ExcessHet.VQSR.vcf.gz",
            apply_sort = false,
            source_vcfs_have_same_header = true,
            create_index = true,
            docker_image = concat_vcfs_docker_image,
            threads = concat_vcfs_threads,
            memory_gb = concat_vcfs_memory_gb
        }

        call step2100_combine_bcftools_stats as step2103_combine_bcftools_stats_all { input:
            sources = step2101_bcftools_stats_all.stats,
            result__name = sub(basename(source_vcf.vcf), ".vcf(.gz)?", "") + ".ExcessHet.VQSR.vcf.gz.bcftools.stats.ALL",
            docker_image = combine_bcftools_stats_docker_image,
            threads = combine_bcftools_stats_threads,
            memory_gb = combine_bcftools_stats_memory_gb
        }

        call step2100_combine_bcftools_stats as step2103_combine_bcftools_stats_pass { input:
            sources = step2102_bcftools_stats_pass.stats,
            result__name = sub(basename(source_vcf.vcf), ".vcf(.gz)?", "") + ".ExcessHet.VQSR.vcf.gz.bcftools.stats.PASS",
            docker_image = combine_bcftools_stats_docker_image,
            threads = combine_bcftools_stats_threads,
            memory_gb = combine_bcftools_stats_memory_gb
        }

    }

    # --------------------------------------------------------------------------------
    # reporting
    # --------------------------------------------------------------------------------

    call multiqc.report as step9001_multiqc_all { input:
        sources = step2103_combine_bcftools_stats_all.result,
        output_prefix = "${analysis_id}.multiqc.ALL",
        docker_image = multiqc_docker_image,
        threads = multiqc_threads,
        memory_gb = multiqc_memory_gb
    }

    call multiqc.report as step9002_multiqc_pass { input:
        sources = step2103_combine_bcftools_stats_pass.result,
        output_prefix = "${analysis_id}.multiqc.PASS",
        docker_image = multiqc_docker_image,
        threads = multiqc_threads,
        memory_gb = multiqc_memory_gb
    }

    call md5sum.md5sum as step9999_md5sum { input:
        sources = flatten([
            [
                step1001_variant_recalibrator_indel.recal_file,
                step1001_variant_recalibrator_indel.recal_index_file,
                step1001_variant_recalibrator_indel.tranches_file,
                step1002_variant_recalibrator_snv.recal_file,
                step1002_variant_recalibrator_snv.recal_index_file,
                step1002_variant_recalibrator_snv.tranches_file,
                step9001_multiqc_all.html,
                step9001_multiqc_all.zip,
                step9002_multiqc_pass.html,
                step9002_multiqc_pass.zip
            ],
            step2004_concat_vcfs.result_vcf_gz,
            step2004_concat_vcfs.result_vcf_gz_index,
            step2103_combine_bcftools_stats_all.result,
            step2103_combine_bcftools_stats_pass.result
        ]),
        md5sum_txt__name = "${analysis_id}.GPCReseq38_0013_Joint_11_VQSR.md5sum.txt",
        docker_image = python_docker_image
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File indel_recal_file = step1001_variant_recalibrator_indel.recal_file
        File indel_recal_index_file = step1001_variant_recalibrator_indel.recal_index_file
        File indel_tranches_file = step1001_variant_recalibrator_indel.tranches_file
        File snv_recal_file = step1002_variant_recalibrator_snv.recal_file
        File snv_recal_index_file = step1002_variant_recalibrator_snv.recal_index_file
        File snv_tranches_file = step1002_variant_recalibrator_snv.tranches_file

        Array[File] vcf_gz = step2004_concat_vcfs.result_vcf_gz
        Array[File] vcf_gz_index = step2004_concat_vcfs.result_vcf_gz_index
        Array[File] bcftools_stats_all = step2103_combine_bcftools_stats_all.result
        Array[File] bcftools_stats_pass = step2103_combine_bcftools_stats_pass.result

        File multiqc_all_html = step9001_multiqc_all.html
        File multiqc_all_zip = step9001_multiqc_all.zip
        File multiqc_pass_html = step9002_multiqc_pass.html
        File multiqc_pass_zip = step9002_multiqc_pass.zip

        File md5sum_txt = step9999_md5sum.md5sum_txt
    }

}


struct ResourceVCF {
    String id
    Boolean known
    Boolean training
    Boolean truth
    Int prior
    File vcf
    File vcf_index
}


struct InputVCF {
    File vcf
    File vcf_index
    String contig
    Int start
    Int end
}


task step0001_filter_by_excess_het {

    input {
        File source_vcf
        File source_vcf_index
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
            --regions ~{region} \
            --targets ~{region} \
            --drop-genotypes \
            --output-type u \
            ~{source_vcf} \
        | bcftools filter \
            --no-version \
            --threads ~{threads} \
            --exclude 'ExcessHet > 54.69' \
            --soft-filter ExcessHet \
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


task step1000_variant_recalibrator {

    input {
        Array[File] vcfs
        Array[File] vcf_indexes
        Array[String] tranches
        Array[String] annotations
        String mode
        Int max_gaussians
        Array[ResourceVCF] resources
        String output_prefix

        String docker_image
        Float memory_gb
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "${memory_gb}G"
    }

    command <<<
        #
        set -euxo pipefail
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        #
        python3.6 <<CODE >resources.txt
        if True:
            import csv
            import sys

            template = '--resource:{id},known={known},training={training},truth={truth},prior={prior} {vcf}'
            with open('~{write_objects(resources)}') as fin:
                sys.stdout.write(' '.join(template.format(**r) for r in csv.DictReader(fin, delimiter='\t')))
        CODE

        #
        gatk VariantRecalibrator \
            ~{sep=" " prefix("--variant ", vcfs)} \
            --trust-all-polymorphic \
            ~{sep=" " prefix("--truth-sensitivity-tranche ", tranches)} \
            ~{sep=" " prefix("--use-annotation ", annotations)} \
            --mode ~{mode} \
            --max-gaussians ~{max_gaussians} \
            $(cat resources.txt) \
            --output ~{output_prefix}.recal \
            --tranches-file ~{output_prefix}.tranches

        #
        rm resources.txt
    >>>

    output {
        File recal_file = "${output_prefix}.recal"
        File recal_index_file = "${output_prefix}.recal.idx"
        File tranches_file = "${output_prefix}.tranches"
    }

}


task step2000_apply_vqsr {

    input {
        File source_vcf
        File source_vcf_index
        File recal_file
        File recal_index_file
        File tranches_file
        Float filter_level
        String mode
        String result_vcf_gz__name

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
        export JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:-} -Xmx~{floor(memory_gb * 1000 * 8/10)}m"

        gatk ApplyVQSR \
            --variant ~{source_vcf} \
            --recal-file ~{recal_file} \
            --tranches-file ~{tranches_file} \
            --truth-sensitivity-filter-level ~{filter_level} \
            --mode ~{mode} \
            --create-output-variant-index true \
            --output ~{result_vcf_gz__name}
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}


task step2003_transfer_filtering_result {

    input {
        File original_vcf
        File original_vcf_index
        File annotation_vcf
        File annotation_vcf_index
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
            --regions ~{region} \
            --targets ~{region} \
            --output-type b \
            --output-file ~{basename(original_vcf)}.temp.bcf.gz \
            ~{original_vcf}

        bcftools index \
            --threads ~{threads} \
            --csi \
            ~{basename(original_vcf)}.temp.bcf.gz

        bcftools annotate \
            --no-version \
            --threads ~{threads} \
            --annotations ~{annotation_vcf} \
            --columns FILTER,INFO/POSITIVE_TRAIN_SITE,INFO/NEGATIVE_TRAIN_SITE,INFO/VQSLOD,INFO/culprit \
            --output-type z \
            --output ~{result_vcf_gz__name} \
            ~{basename(original_vcf)}.temp.bcf.gz

        bcftools index \
            --threads ~{threads} \
            --tbi \
            ~{result_vcf_gz__name}

        rm ~{basename(original_vcf)}.temp.bcf.gz*
    >>>

    output {
        File result_vcf_gz = "${result_vcf_gz__name}"
        File result_vcf_gz_index = "${result_vcf_gz__name}.tbi"
    }

}


task step2100_bcftools_stats {

    input {
        File source_vcf
        File source_vcf_index
        String? filter
        String stats__name

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

        bcftools stats \
            --threads ~{threads} \
            --samples-file <(bcftools query --list-samples ~{source_vcf}) \
            ~{"--apply-filters " + filter} \
            ~{source_vcf} \
        > ~{stats__name}
    >>>

    output {
        File stats = "${stats__name}"
    }

}


task step2100_combine_bcftools_stats {

    input {
        Array[File] sources
        String result__name

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

        python3 - <<EOS >~{result__name}
        if True:
            def _read_bcftools_stats(fin):
                sn = {}
                tstv = None
                sis = None
                idd = {}
                st = {}
                psc = {}
                psi = {}

                for line in fin:
                    #
                    line = line.strip()
                    if (not line) or line.startswith('#'):
                        continue

                    #
                    cols = line.split('\t')
                    assert cols[1] == '0'

                    if cols[0] in ('ID', 'AF', 'QUAL', 'DP', 'HWE'):
                        pass    # unsupported
                    elif line.startswith('SN'):
                        sn[cols[2][:-1]] = int(cols[3])
                    elif line.startswith('TSTV'):
                        tstv = tuple(t(c) for t, c in zip((int, int, float, int, int, float), cols[2:8]))
                    elif line.startswith('SiS'):
                        sis = tuple(map(int, cols[2:]))
                    elif line.startswith('IDD'):
                        idd[cols[2]] = int(cols[3])
                    elif line.startswith('ST'):
                        st[cols[2]] = int(cols[3])
                    elif line.startswith('PSC'):
                        psc[cols[2]] = tuple(t(c) for t, c in zip((int, int, int, int, int, int, float, int, int, int), cols[3:13]))
                    elif line.startswith('PSI'):
                        psi[cols[2]] = tuple(t(c) for t, c in zip((int, int, int, float, int, int, int), cols[3:9]))
                    else:
                        raise Exception(line.split('\t')[0])

                return sn, tstv, sis, idd, st, psc, psi

            def _merge_bcftools_stats(records):
                #
                if len(set(r[0]['number of samples'] for r in records)) != 1:
                    raise Exception

                #
                sn_all = {k: 0 for k in records[0][0]}
                tstv_all = 0, 0, 0.0, 0, 0, 0.0
                sis_all = 0, 0, 0, 0, 0, 0, 0, 0
                idd_all = {}
                st_all = {}
                psc_all = {s: (0, 0, 0, 0, 0, 0, 0.0, 0, 0, 0) for s in records[0][5]}
                psi_all = {s: (0, 0, 0, 0.0, 0, 0) for s in records[0][6]}

                for sn, tstv, sis, idd, st, psc, psi in records:
                    sn_all = {k: (sn_all[k] + sn[k]) for k in sn}
                    tstv_all = [x + y for x, y in zip(tstv_all, tstv)]
                    sis_all = [x + y for x, y in zip(sis_all, sis)]
                    psc_all = {s: [x + y for x, y in zip(psc_all[s], psc[s])] for s in psc}
                    psi_all = {s: [x + y for x, y in zip(psi_all[s], psi[s])] for s in psi}
                    for key in idd:
                        idd_all[key] = idd_all.get(key, 0) + idd[key]
                    for key in st:
                        st_all[key] = st_all.get(key, 0) + st[key]

                sn_all['number of samples'] = records[0][0]['number of samples']
                tstv_all[2] = tstv_all[0] / tstv_all[1]     # ts/tv
                tstv_all[5] = tstv_all[3] / tstv_all[4]     # ts/tv
                for sample in psc_all:
                    psc_all[sample][6] /= len(records)      # average depth
                for sample in psi_all:
                    try:
                        psi_all[sample][3] = psi_all[sample][1] / (psi_all[sample][0] + psi_all[sample][1]) # out/(in+out) ratio
                    except ZeroDivisionError:
                        psi_all[sample][3] = 0.0

                return sn_all, tstv_all, sis_all, idd_all, st_all, psc_all, psi_all

            def _main():
                records = []
                for line in open('~{write_lines(sources)}'):
                    with open(line.strip()) as fin:
                        records.append(_read_bcftools_stats(fin))

                sn, tstv, sis, idd, st, psc, psi = _merge_bcftools_stats(records)

                print('# This file was produced by bcftools stats and can be plotted using plot-vcfstats.')
                print('\t'.join(['ID', '0', '~{result__name}']))

                print('# SN\t[2]id\t[3]key\t[4]value')
                for key, value in sn.items():
                    print(f'SN\t0\t{key}\t{value}')

                print('# SiS, Singleton stats:')
                print('\t'.join([
                    '# SiS', '[2]id', '[3]allele count', '[4]number of SNPs',
                    '[5]number of transitions', '[6]number of transversions',
                    '[7]number of indels', '[8]repeat-consistent', '[9]repeat-inconsistent', '[10]not applicable'
                ]))
                print('\t'.join(['SiS', '0'] + list(map(str, sis))))

                print('# TSTV, transitions/transversions:')
                print('\t'.join(['# TSTV', '[2]id', '[3]ts', '[4]tv', '[5]ts/tv', '[6]ts (1st ALT)', '[7]tv (1st ALT)', '[8]ts/tv (1st ALT)']))
                print('\t'.join(['TSTV', '0'] + list(map(str, tstv))))

                print('# SiS, transitions/transversions:')
                print('\t'.join([
                    '# SiS', '[2]id', '[3]allele count', '[4]number of SNPs',
                    '[5]number of transitions', '[6]number of transversions',
                    '[7]number of indels', '[8]repeat-consistent', '[9]repeat-inconsistent', '[10]not applicable'
                ]))
                print('\t'.join(['SiS', '0'] + list(map(str, sis))))

                print('# IDD, InDel distribution:')
                print('\t'.join(['# IDD', '[2] id', '[3]length (deletions negative)', '[4]count']))
                for length, count in idd.items():
                    print('\t'.join(['IDD', '0', length, str(count)]))

                print('# ST, Substitution types:')
                print('\t'.join(['# IDD', '[2] id', '[3]type', '[4]count']))
                for type, count in st.items():
                    print('\t'.join(['ST', '0', type, str(count)]))

                print('# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. Haploid counts include both SNPs and indels.')
                print('\t'.join([
                    '# PSC', '[2]id', '[3]sample', '[4]nRefHom', '[5]nNonRefHom', '[6]nHets',
                    '[7]nTransitions', '[8]nTransversions', '[9]nIndels', '[10]average depth', '[11]nSingletons',
                    '[12]nHapRef', '[13]nHapAlt'
                ]))
                for sample, values in psc.items():
                    print('\t'.join(['PSC', '0', sample] + list(map(str, values))))

                print('# PSI, Per-Sample Indels')
                print('\t'.join(['# PSI', '[2]id', '[3]sample', '[4]in-frame', '[5]out-frame', '[6]not applicable', '[7]out/(in+out) ratio', '[8]nHets', '[9]nAA']))
                for sample, values in psi.items():
                    print('\t'.join(['PSI', '0', sample] + list(map(str, values))))

            #
            _main()
        EOS
    >>>

    output {
        File result = "${result__name}"
    }

}
