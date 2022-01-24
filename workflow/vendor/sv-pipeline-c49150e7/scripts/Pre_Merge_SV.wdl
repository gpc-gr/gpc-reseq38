version 1.0
import "Pre_Merge_SV_per_sample.wdl" as per_sample
import "Pre_Merge_QC_per_sample.wdl" as qc
import "SV_Tasks.wdl" as SV

workflow Pre_Merge_SV {
  input {
    Array[File] aligned_crams
    String aligned_cram_suffix

    # reference inputs
    File ref_fasta
    File ref_fasta_index
    File ref_cache
    File? call_regions_bed
    File? call_regions_bed_index
    File exclude_regions
    String cohort
    String center

    # system inputs
    Int preemptible_tries
  }


  scatter (i in range(length(aligned_crams))) {
    File aligned_cram = aligned_crams[i]

    call per_sample.Pre_Merge_SV_Per_Sample {
      input:
        aligned_cram = aligned_cram,
        aligned_cram_suffix = aligned_cram_suffix,
	    ref_fasta = ref_fasta,
	    ref_fasta_index = ref_fasta_index,
        call_regions_bed = call_regions_bed,
        call_regions_bed_index = call_regions_bed_index,
	    ref_cache = ref_cache,
	    exclude_regions = exclude_regions,
	    preemptible_tries = preemptible_tries
    }
    
    call qc.Pre_Merge_QC_Per_Sample {
      input:
        manta_vcf = Pre_Merge_SV_Per_Sample.manta_vcf,
        lumpy_vcf = Pre_Merge_SV_Per_Sample.smoove_vcf,
        cnvnator_vcf = Pre_Merge_SV_Per_Sample.cnvnator_output_cn_txt,
        cohort = cohort,
        center = center,
	preemptible_tries = preemptible_tries
    }
  }

  #scatter (p in [("manta", Pre_Merge_QC_Per_Sample.manta_counts), ("lumpy", Pre_Merge_QC_Per_Sample.lumpy_counts)]) {
  #  call SV.Make_Count_Plot {
  #    input:
  #      name=p.left,
  #      count_files=p.right
  #  }
  #}

  output {
    Array[File] cram_indices = Pre_Merge_SV_Per_Sample.cram_index
    Array[File] manta_vcfs = Pre_Merge_SV_Per_Sample.manta_vcf
    Array[File] manta_tbis = Pre_Merge_SV_Per_Sample.manta_tbi
    Array[File] manta_original_vcfs = Pre_Merge_SV_Per_Sample.manta_original_vcf
    Array[File] manta_original_tbis = Pre_Merge_SV_Per_Sample.manta_original_tbi
    Array[File] cnvnator_cn_hist_roots = Pre_Merge_SV_Per_Sample.cnvnator_cn_hist_root
    Array[File] cnvnator_output_cn_txt_files = Pre_Merge_SV_Per_Sample.cnvnator_output_cn_txt
    Array[File] cnvnator_cn_bed_files = Pre_Merge_SV_Per_Sample.cnvnator_cn_bed
    Array[File] smoove_vcfs = Pre_Merge_SV_Per_Sample.smoove_vcf
    Array[File] smoove_csis = Pre_Merge_SV_Per_Sample.smoove_csi
    Array[File] lumpy_counts = Pre_Merge_QC_Per_Sample.lumpy_counts
    Array[File] manta_counts = Pre_Merge_QC_Per_Sample.manta_counts
    #Array[File] count_plots = Make_Count_Plot.counts_plot
  }
}
