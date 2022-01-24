#
#
#

version 1.0

import "./modules/vcf2plink.wdl"


workflow GPCReseq38_0013_Joint_19_PCA {

    # --------------------------------------------------------------------------------
    # input
    # --------------------------------------------------------------------------------

    input {
        String analysis_id

        File reference_fasta
        Array[File] reference_fasta_general_indexes = [
            "${reference_fasta}.fai",
            sub(reference_fasta, ".fa(sta)?$", ".dict")
        ]

        Array[Dataset] datasets

        Int chunk_size = 3 * 1000 * 1000    # 3Mb

        String bcftools_docker_image = "quay.io/biocontainers/bcftools:1.11--h7c999a4_0"
        String flashpca_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-flashpca:flashpca_2.0"
        String plink_docker_image = "quay.io/biocontainers/plink:1.90b6.21--h779adbc_1"
        String python_docker_image = "python:3.8.6-slim-buster"
        String r_docker_image = "registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-r:dplyr_1.0.7-ggplot2_3.3.5"

        String split_region_docker_image = python_docker_image
        Float split_region_memory_gb = 1
        String make_slim_vcf_docker_image = bcftools_docker_image
        Int make_slim_vcf_threads = 4
        Float make_slim_vcf_memory_gb = 16
        String concat_vcfs_docker_image = bcftools_docker_image
        Int concat_vcfs_threads = 8
        Float concat_vcfs_memory_gb = 16
        String vcf2plink_docker_image = plink_docker_image
        Int vcf2plink_threads = 1
        Float vcf2plink_memory_gb = 90
        String update_sex_in_fam_docker_image = python_docker_image
        Float update_sex_in_fam_memory_gb = 4

        String collect_common_variants_docker_image = python_docker_image
        Int collect_common_variants_threads = 1
        Float collect_common_variants_memory_gb = 4
        String extract_samples_and_variants_docker_image = plink_docker_image
        Int extract_samples_and_variants_threads = 1
        Float extract_samples_and_variants_memory_gb = 32
        String filter_idtable_docker_image = python_docker_image
        Int filter_idtable_threads = 1
        Float filter_idtable_memory_gb = 4
        String combine_datasets_docker_image = plink_docker_image
        Int combine_datasets_threads = 1
        Float combine_datasets_memory_gb = 90
        String combine_idtables_docker_image = python_docker_image
        Int combine_idtables_threads = 4
        Float combine_idtables_memory_gb = 4
        String combined_qc_docker_image = plink_docker_image
        Int combined_qc_threads = 1
        Float combined_qc_memory_gb = 32

        String pruning_docker_image = plink_docker_image
        Int pruning_threads = 8
        Float pruning_memory_gb = 32
        String pca_docker_image = flashpca_docker_image
        Int pca_threads = 16
        Float pca_memory_gb = 64
        String plot_eigenvec_docker_image = r_docker_image
        Int plot_eigenvec_threads = 1
        Float plot_eigenvec_memory_gb = 4
    }

    # --------------------------------------------------------------------------------
    # combine datasets
    # --------------------------------------------------------------------------------

    scatter (dataset in datasets) {

        scatter (entry in dataset.vcfs) {
            File source_vcfs = entry.vcf
            File source_vcf_indexes = entry.vcf_index
            String contigs = entry.contig
            Int starts = entry.start
            Int ends = entry.end
        }

        call vcf2plink.vcf2plink as step0001_vcf2plink { input:
            reference_fasta = reference_fasta,
            reference_fasta_general_indexes = reference_fasta_general_indexes,
            source_vcfs = source_vcfs,
            source_vcf_indexes = source_vcf_indexes,
            contigs = contigs,
            starts = starts,
            ends = ends,
            output_prefix = "${analysis_id}.QCed",
            target_sample_list = dataset.target_sample_list,
            pass_only = true,
            snv_only = true,
            biallelic_only = true,
            hwe_ge = 0.0001,
            maf_ge = 0.01,
            geno_le = 0.01,
            normalize_position_and_alleles = true,
            chunk_size = chunk_size,
            sample_idtable = dataset.sample_idtable,
            sample_sex_column_key = dataset.sample_sex_column_key,
            split_region_docker_image = split_region_docker_image,
            split_region_memory_gb = split_region_memory_gb,
            make_slim_vcf_docker_image = make_slim_vcf_docker_image,
            make_slim_vcf_threads = make_slim_vcf_threads,
            make_slim_vcf_memory_gb = make_slim_vcf_memory_gb,
            concat_vcfs_docker_image = concat_vcfs_docker_image,
            concat_vcfs_threads = concat_vcfs_threads,
            concat_vcfs_memory_gb = concat_vcfs_memory_gb,
            vcf2plink_docker_image = vcf2plink_docker_image,
            vcf2plink_threads = vcf2plink_threads,
            vcf2plink_memory_gb = vcf2plink_memory_gb,
            update_sex_in_fam_docker_image = update_sex_in_fam_docker_image,
            update_sex_in_fam_memory_gb = update_sex_in_fam_memory_gb
        }

    }

    call step1001_collect_common_variants { input:
        bims = step0001_vcf2plink.bim,
        result_txt__name = "${analysis_id}.common_variants.txt",
        docker_image = collect_common_variants_docker_image,
        threads = collect_common_variants_threads,
        memory_gb = collect_common_variants_memory_gb
    }

    scatter (index in range(length(datasets))) {

        call step1002_extract_variants { input:
            source_bed = step0001_vcf2plink.bed[index],
            source_bim = step0001_vcf2plink.bim[index],
            source_fam = step0001_vcf2plink.fam[index],
            target_variants_txt = step1001_collect_common_variants.result_txt,
            output_prefix = sub(basename(step0001_vcf2plink.bed[index]), ".bed$", ".cleaned.extracted"),
            docker_image = extract_samples_and_variants_docker_image,
            threads = extract_samples_and_variants_threads,
            memory_gb = extract_samples_and_variants_memory_gb
        }

        call step1003_filter_idtable { input:
            source_idtable = datasets[index].sample_idtable,
            target_sample_list = datasets[index].target_sample_list,
            sample_category_column_key = datasets[index].sample_category_column_key,
            result_idtable__name = sub(basename(datasets[index].sample_idtable), ".tsv$", ".extracted.tsv"),
            docker_image = filter_idtable_docker_image,
            threads = filter_idtable_threads,
            memory_gb = filter_idtable_memory_gb
        }

    }

    call step1008_combine_datasets { input:
        source_beds = step1002_extract_variants.result_bed,
        source_bims = step1002_extract_variants.result_bim,
        source_fams = step1002_extract_variants.result_fam,
        output_prefix = "${analysis_id}.combined.genotype.raw",
        docker_image = combine_datasets_docker_image,
        threads = combine_datasets_threads,
        memory_gb = combine_datasets_memory_gb
    }

    call step1009_combine_idtables { input:
        source_idtables = step1003_filter_idtable.result_idtable,
        result_idtable__name = "${analysis_id}.combined.itable.tsv",
        docker_image = combine_idtables_docker_image,
        threads = combine_idtables_threads,
        memory_gb = combine_idtables_memory_gb
    }

    # --------------------------------------------------------------------------------
    # PCA
    # --------------------------------------------------------------------------------

    call step2001_combined_qc { input:
        source_bed = step1008_combine_datasets.result_bed,
        source_bim = step1008_combine_datasets.result_bim,
        source_fam = step1008_combine_datasets.result_fam,
        output_prefix = "${analysis_id}.combined.genotype.clenaed",
        docker_image = combined_qc_docker_image,
        threads = combined_qc_threads,
        memory_gb = combined_qc_memory_gb
    }

    call step2002_pruning { input:
        source_bed = step2001_combined_qc.result_bed,
        source_bim = step2001_combined_qc.result_bim,
        source_fam = step2001_combined_qc.result_fam,
        result_prune_file_prefix = "${analysis_id}.combined.genotype.cleaned",
        result_bfile_prefix = "${analysis_id}.combined.genotype.cleaned.pruned",
        docker_image = pruning_docker_image,
        threads = pruning_threads,
        memory_gb = pruning_memory_gb
    }

    call step2003_pca { input:
        bed = step2002_pruning.pruned_bed,
        bim = step2002_pruning.pruned_bim,
        fam = step2002_pruning.pruned_fam,
        output_prefix = "${analysis_id}.combined.genotype.cleaned.pruned.pca",
        docker_image = pca_docker_image,
        threads = pca_threads,
        memory_gb = pca_memory_gb
    }

    call step2004_plot_eigenvec { input:
        eigenvec_txt = step2003_pca.eigenvec_txt,
        idtable = step1009_combine_idtables.result_idtable,
        eigenvec_pdf__name = "${basename(step2003_pca.eigenvec_txt)}.pdf",
        docker_image = plot_eigenvec_docker_image,
        threads = plot_eigenvec_threads,
        memory_gb = plot_eigenvec_memory_gb
    }

    # --------------------------------------------------------------------------------
    # output
    # --------------------------------------------------------------------------------

    output {
        File combined_idtable = step1009_combine_idtables.result_idtable

        Array[File] raw_dataset = [
            step1008_combine_datasets.result_bed,
            step1008_combine_datasets.result_bim,
            step1008_combine_datasets.result_fam,
        ]
        Array[File] cleaned_dataset = [
            step2001_combined_qc.result_bed,
            step2001_combined_qc.result_bim,
            step2001_combined_qc.result_fam
        ]
        Array[File] pruned_dataset = [
            step2002_pruning.prune_in,
            step2002_pruning.prune_out,
            step2002_pruning.pruned_bed,
            step2002_pruning.pruned_bim,
            step2002_pruning.pruned_fam
        ]

        Array[File] pca_results = [
            step2003_pca.eigenvec_txt,
            step2003_pca.pcs_txt,
            step2003_pca.eigenval_txt,
            step2003_pca.pve_txt,
            step2004_plot_eigenvec.eigenvec_pdf
        ]
    }

}


struct VCF {
    File vcf
    File vcf_index
    String contig
    Int start
    Int end
}


struct Dataset {
    Array[VCF] vcfs
    File sample_idtable
    String sample_sex_column_key
    String sample_category_column_key
    File target_sample_list
}


task step0001_individual_qc {

    input {
        File source_bed
        File source_bim
        File source_fam
        File target_samples_txt
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
        set -euxo pipefail

        cat ~{source_fam} | grep -w -f ~{target_samples_txt} > keep.fam

        plink \
            --bed ~{source_bed} \
            --bim ~{source_bim} \
            --fam ~{source_fam} \
            --keep keep.fam \
            --maf 0.01 \
            --hwe 1e-5 \
            --make-just-bim \
            --out ~{output_prefix} \
            --memory ~{floor(memory_gb * 1000 * 8/10)}

        rm ~{output_prefix}.temp1.* keep.fam
    >>>

    output {
        File result_bim = "${output_prefix}.bim"
    }

}


task step1001_collect_common_variants {

    input {
        Array[File] bims
        String result_txt__name

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

        python3 - <<EOS > ~{result_txt__name}
        if True:
            import collections

            variants = collections.defaultdict(int)
            for source in open('~{write_lines(bims)}'):
                with open(source) as fin:
                    for line in fin:
                        id = line.strip().split()[1]
                        variants[id] += 1

            for key, count in variants.items():
                if count == ~{length(bims)}:
                    print(key)
        EOS
    >>>

    output {
        File result_txt = "${result_txt__name}"
    }

}


task step1002_extract_variants {

    input {
        File source_bed
        File source_bim
        File source_fam
        File target_variants_txt
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
        set -euxo pipefail

        plink \
            --bed ~{source_bed} \
            --bim ~{source_bim} \
            --fam ~{source_fam} \
            --extract ~{target_variants_txt}
            --make-bed \
            --out ~{output_prefix} \
            --memory ~{floor(memory_gb * 1000 * 8/10)}

        rm keep.fam
    >>>

    output {
        File result_bed = "${output_prefix}.bed"
        File result_bim = "${output_prefix}.bim"
        File result_fam = "${output_prefix}.fam"
    }

}


task step1003_filter_idtable {

    input {
        File source_idtable
        File target_sample_list
        String sample_category_column_key
        String result_idtable__name

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
        set -euxo pipeline

        python3 - <<EOS > ~{result_idtable__name}
        if True:
            import csv

            with open('~{target_sample_list}') as fin:
                targets = set(l.strip() for l in fin if l.strip())

            with open('~{source_idtable}') as fin:
                print('id\tcategory')
                for record in csv.DictReader(fin, delimiter='\t'):
                    if record['id'] in targets:
                        print('{}\t{}'.format(record['id'], record['~{sample_category_column_key}']))
        EOS
    >>>

    output {
        File result_idtable = "${result_idtable__name}"
    }

}


task step1008_combine_datasets {

    input {
        Array[File] source_beds
        Array[File] source_bims
        Array[File] source_fams
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
        set -euxo pipefail

        paste ~{write_lines(source_beds)} ~{write_lines(source_bims)} ~{write_lines(source_fams)} > sources.txt

        plink \
            --merge-list sources.txt \
            --make-bed \
            --out ~{output_prefix} \
            --memory ~{floor(memory_gb * 1000 * 8/10)}
    >>>

    output {
        File result_bed = "${output_prefix}.bed"
        File result_bim = "${output_prefix}.bim"
        File result_fam = "${output_prefix}.fam"
    }

}


task step1009_combine_idtables {

    input {
        Array[File] source_idtables
        String result_idtable__name

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

        {
            head -n1 ~{source_idtables[0]}
            echo ~{sep=' ' source_idtables} | xargs -n1 awk 'NR > 1'
        } > ~{result_idtable__name}
    >>>

    output {
        File result_idtable = "${result_idtable__name}"
    }

}


task step2001_combined_qc {

    input {
        File source_bed
        File source_bim
        File source_fam
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
        set -euxo pipefail

        plink \
            --bed ~{source_bed} \
            --bim ~{source_bim} \
            --fam ~{source_fam} \
            --maf 0.05 \
            --hwe 0.05 \
            --geno 0.01 \
            --make-bed \
            --out ~{output_prefix} \
            --memory ~{floor(memory_gb * 1000 * 8/10)}
    >>>

    output {
        File result_bed = "${output_prefix}.bed"
        File result_bim = "${output_prefix}.bim"
        File result_fam = "${output_prefix}.fam"
    }

}


task step2002_pruning {

    input {
        File source_bed
        File source_bim
        File source_fam
        String result_prune_file_prefix
        String result_bfile_prefix

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

        plink \
            --bed ~{source_bed} \
            --bim ~{source_bim} \
            --fam ~{source_fam} \
            --indep-pairwise 1500 150 0.03 \
            --out ~{result_prune_file_prefix} \
            --threads ~{threads} \
            --memory ~{floor(memory_gb * 1000 * 8/10)}

        plink \
            --bed ~{source_bed} \
            --bim ~{source_bim} \
            --fam ~{source_fam} \
            --extract ~{result_prune_file_prefix}.prune.in \
            --out ~{result_bfile_prefix} \
            --memory ~{floor(memory_gb * 1000 * 8/10)}
    >>>

    output {
        File prune_in = "${result_prune_file_prefix}.prune.in"
        File prune_out = "${result_prune_file_prefix}.prune.out"
        File pruned_bed = "${result_bfile_prefix}.bed"
        File pruned_bim = "${result_bfile_prefix}.bim"
        File pruned_fam = "${result_bfile_prefix}.fam"
    }

}


task step2003_pca {

    input {
        File bed
        File bim
        File fam
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
        set -euxo pipefail

        flashpca \
            --verbose \
            --seed 1 \
            --bed ~{bed} \
            --bim ~{bim} \
            --fam ~{fam} \
            --numthreads ~{threads}

        mv eigenvectors.txt ~{output_prefix}.eigenvectors.txt
        mv pcs.txt ~{output_prefix}.pcs.txt
        mv eigenvalues.txt ~{output_prefix}.eigenvalues.txt
        mv pve.txt ~{output_prefix}.pve.txt
    >>>

    output {
        File eigenvec_txt = "${output_prefix}.eigenvectors.txt"
        File pcs_txt = "${output_prefix}.pcs.txt"
        File eigenval_txt = "${output_prefix}.eigenval.txt"
        File pve_txt = "${output_prefix}.pve.txt"
    }

}


task step2004_plot_eigenvec {

    input {
        File eigenvec_txt
        File idtable
        String eigenvec_pdf__name

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

        Rscript - <<EOS
            library(dplyr)
            library(ggplot2)

            d.eigenvec <- read.table('~{eigenvec_txt}', head=T)
            d.idtable <- read.table('~{idtable}', head=T)
            d <- left_join(d.eigenvec, d.idtable, by=c(id='id'))

            pdf('~{eigenvec_pdf__name}')
            plot(ggplot(d) +
                geom_point(aes(x=PC1, y=PC2, colour=category)) +
                ggtitle('~{basename(eigenvec_txt)} :: PC1 vs PC2') +
                xlab('PC1') +
                xlab('PC2')
            )
            plot(ggplot(d) +
                geom_point(aes(x=PC1, y=PC2, colour=category)) +
                ggtitle('~{basename(eigenvec_txt)} :: PC1 vs PC2') +
                xlab('PC1') +
                xlab('PC3')
            )
        EOS
    >>>

    output {
        File eigenvec_pdf = "${eigenvec_pdf__name}"
    }

}
