version 1.0

import "GPCReseq38_0011_SingleSample_01_Align.wdl" as Align
import "GPCReseq38_0011_SingleSample_02_BQSR.wdl" as BQSR
import "GPCReseq38_0011_SingleSample_11_Autosome.wdl" as Autosome
import "GPCReseq38_0011_SingleSample_12_chrXY_PAR2.wdl" as chrXY_PAR2
import "GPCReseq38_0011_SingleSample_13_chrXY_PAR3.wdl" as chrXY_PAR3
import "GPCReseq38_0011_SingleSample_14_Mitochondria.wdl" as Mitochondria

workflow GPCReseq38_0011_SingleSample_all {
    input {
        # Initial alignment
        File reference_fasta
        Array[File] reference_fasta_general_indexes = [
        "${reference_fasta}.fai",
        sub(reference_fasta, ".fa(sta)?$", ".dict")
        ]
        Array[File] reference_fasta_bwa_indexes = [
        "${reference_fasta}.64.amb",
        "${reference_fasta}.64.ann",
        "${reference_fasta}.64.bwt",
        "${reference_fasta}.64.pac",
        "${reference_fasta}.64.sa",
        "${reference_fasta}.64.alt"
        ]

        Array[String] wgs_metrics_regions = [
        "autosome",
        "chrXY_PAR2.chrX_nonPAR",
        "chrXY_PAR2.chrY_nonPAR",
        "mitochondria"
        ]

        String sample_id
        Array[ReadPair] sample_read_pairs

        String? accelaration

        # BQSR
        String bqsr_model_region_type = "autosome"
        File bqsr_model_region_list = "${reference_fasta}.regions.${bqsr_model_region_type}.tsv"

        Array[VCF] known_site_vcfs

        # PAR3
        File chrXY_PAR3_reference_fasta
        Array[File] chrXY_PAR3_reference_fasta_general_indexes = [
            "${chrXY_PAR3_reference_fasta}.fai",
            sub(chrXY_PAR3_reference_fasta, ".fa(sta)?$", ".dict")
        ]
        Array[File] chrXY_PAR3_reference_fasta_bwa_indexes = [
            "${chrXY_PAR3_reference_fasta}.64.amb",
            "${chrXY_PAR3_reference_fasta}.64.ann",
            "${chrXY_PAR3_reference_fasta}.64.bwt",
            "${chrXY_PAR3_reference_fasta}.64.pac",
            "${chrXY_PAR3_reference_fasta}.64.sa",
            "${chrXY_PAR3_reference_fasta}.64.alt"
        ]

        # Mitocondria
        File mitochondria_shifted_reference_fasta
        Array[File] mitochondria_shifted_reference_fasta_general_indexes = [
            "${mitochondria_shifted_reference_fasta}.fai",
            sub(mitochondria_shifted_reference_fasta, ".fa(sta)?$", ".dict")
        ]
        Array[File] mitochondria_shifted_reference_fasta_bwa_indexes = [
            "${mitochondria_shifted_reference_fasta}.64.amb",
            "${mitochondria_shifted_reference_fasta}.64.ann",
            "${mitochondria_shifted_reference_fasta}.64.bwt",
            "${mitochondria_shifted_reference_fasta}.64.pac",
            "${mitochondria_shifted_reference_fasta}.64.sa",
            "${mitochondria_shifted_reference_fasta}.64.alt"
        ]
    }

    call Align.GPCReseq38_0011_SingleSample_01_Align {
        input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        reference_fasta_bwa_indexes = reference_fasta_bwa_indexes,
        wgs_metrics_regions = wgs_metrics_regions,
        sample_id = sample_id,
        sample_read_pairs = sample_read_pairs,
        accelaration = accelaration
    }

    call BQSR.GPCReseq38_0011_SingleSample_02_BQSR {
        input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        bqsr_model_region_type = bqsr_model_region_type,
        bqsr_model_region_list = bqsr_model_region_list,
        known_site_vcfs = known_site_vcfs,
        sample_id = sample_id,
        sample_raw_cram = GPCReseq38_0011_SingleSample_01_Align.cram,
        sample_raw_cram_index = GPCReseq38_0011_SingleSample_01_Align.cram_index
    }

    call Autosome.GPCReseq38_0011_SingleSample_11_Autosome {
        input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        sample_id = sample_id,
        sample_cram = GPCReseq38_0011_SingleSample_01_Align.cram,
        sample_cram_index = GPCReseq38_0011_SingleSample_01_Align.cram_index,
        accelaration = accelaration
    }

    call chrXY_PAR2.GPCReseq38_0011_SingleSample_12_chrXY_PAR2 {
        input:
        reference_fasta = reference_fasta,
        reference_fasta_general_indexes = reference_fasta_general_indexes,
        sample_id = sample_id,
        sample_cram = GPCReseq38_0011_SingleSample_01_Align.cram,
        sample_cram_index = GPCReseq38_0011_SingleSample_01_Align.cram_index,
        accelaration = accelaration
    }

    call chrXY_PAR3.GPCReseq38_0011_SingleSample_13_chrXY_PAR3 {
        input:
        baseline_reference_fasta = reference_fasta,
        baseline_reference_fasta_general_indexes = reference_fasta_general_indexes,
        chrXY_PAR3_reference_fasta = chrXY_PAR3_reference_fasta, 
        chrXY_PAR3_reference_fasta_general_indexes = chrXY_PAR3_reference_fasta_general_indexes,
        chrXY_PAR3_reference_fasta_bwa_indexes = chrXY_PAR3_reference_fasta_bwa_indexes,
        sample_id = sample_id,
        sample_baseline_cram = GPCReseq38_0011_SingleSample_01_Align.cram,
        sample_baseline_cram_index = GPCReseq38_0011_SingleSample_01_Align.cram_index,
        accelaration = accelaration
    }

    call Mitochondria.GPCReseq38_0011_SingleSample_14_Mitochondria {
        input:
        baseline_reference_fasta = reference_fasta,
        baseline_reference_fasta_general_indexes = reference_fasta_general_indexes,
        mitochondria_shifted_reference_fasta = mitochondria_shifted_reference_fasta,
        mitochondria_shifted_reference_fasta_general_indexes = mitochondria_shifted_reference_fasta_general_indexes,
        mitochondria_shifted_reference_fasta_bwa_indexes = mitochondria_shifted_reference_fasta_bwa_indexes,
        sample_id = sample_id,
        sample_baseline_cram = GPCReseq38_0011_SingleSample_01_Align.cram,
        sample_baseline_cram_index = GPCReseq38_0011_SingleSample_01_Align.cram_index,
        accelaration = accelaration
    }


    output {
        # Align
        Array[File] read_metrics = GPCReseq38_0011_SingleSample_01_Align.read_metrics

        File cram = GPCReseq38_0011_SingleSample_01_Align.cram
        File cram_index = GPCReseq38_0011_SingleSample_01_Align.cram_index
        Array[File] cram_metrics = GPCReseq38_0011_SingleSample_01_Align.cram_metrics

        File multiqc_html = GPCReseq38_0011_SingleSample_01_Align.multiqc_html
        File multiqc_zip = GPCReseq38_0011_SingleSample_01_Align.multiqc_zip

        File md5sum_txt_align = GPCReseq38_0011_SingleSample_01_Align.md5sum_txt

        # BQSR
        File bqsr_table = GPCReseq38_0011_SingleSample_02_BQSR.bqsr_table
        File md5sum_txt_bqsr = GPCReseq38_0011_SingleSample_02_BQSR.md5sum_txt

        # Autosome call
        File autosome_gvcf_gz = GPCReseq38_0011_SingleSample_11_Autosome.gvcf_gz
        File autosome_gvcf_gz_index = GPCReseq38_0011_SingleSample_11_Autosome.gvcf_gz_index

        File md5sum_txt_autosome = GPCReseq38_0011_SingleSample_11_Autosome.md5sum_txt

        # chrXY PAR2 call
        Array[File] chrXY_PAR2_gvcf_gz = GPCReseq38_0011_SingleSample_12_chrXY_PAR2.gvcf_gz
        Array[File] chrXY_PAR2_gvcf_gz_index = GPCReseq38_0011_SingleSample_12_chrXY_PAR2.gvcf_gz_index

        File md5sum_txt_chrXY_PAR2 = GPCReseq38_0011_SingleSample_12_chrXY_PAR2.md5sum_txt

        # chrXY PAR3 call
        Array[File] chrXY_PAR3_gvcf_gz = GPCReseq38_0011_SingleSample_13_chrXY_PAR3.gvcf_gz
        Array[File] chrXY_PAR3_gvcf_gz_index = GPCReseq38_0011_SingleSample_13_chrXY_PAR3.gvcf_gz_index

        # Mitochondria.
        File mitochondria_shifted_cram = GPCReseq38_0011_SingleSample_14_Mitochondria.shifted_cram
        File mitochondria_shifted_cram_index = GPCReseq38_0011_SingleSample_14_Mitochondria.shifted_cram_index
        File mitochondria_shifted_gvcf_gz = GPCReseq38_0011_SingleSample_14_Mitochondria.shifted_gvcf_gz
        File mitochondria_shifted_gvcf_gz_index = GPCReseq38_0011_SingleSample_14_Mitochondria.shifted_gvcf_gz_index

    }
}
