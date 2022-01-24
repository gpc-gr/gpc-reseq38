gpc-reseq38
===========

Workflows
---------
| name                                       | summary                                       |
| ------------------------------------------ | --------------------------------------------- |
| 0001_ReferencePreparation                  | Prepares index files                          |
| 0011_SingleSample_all                      | Runs all 0011_SingleSample_\d+ at single run  |
| 0011_SingleSample_01_Align                 | Aligns FASTQs onto a reference genome         |
| 0011_SingleSample_02_BQSR                  | Generates BQSR table                          |
| 0011_SingleSample_11_Autosome              | Call variants on autosome                     |
| 0011_SingleSample_12_chrXY_PAR2            | Call variants on chrX and chrY (PAR2 version) |
| 0011_SingleSample_13_chrXY_PAR3            | Call variants on chrX and chrY (PAR3 version) |
| 0011_SingleSample_14_Mitochondria          | Call variants on mitochondria                 |
| 0011_SingleSample_21_SVPreMerge            | Runs sv-pipeline (Pre_Merge_SV_per_sample)    |
| 0011_SingleSample_91_VerifyBAMID           | Runs verifyBAMID2 to check sample swap        |
| 0011_SingleSample_92_DepthDetail           | Collects depth information from CRAM          |
| 0011_SingleSample_95_GVCF2VCF_Autosome     | Runs GenotypeGVCFs (autosome)                 |
| 0011_SingleSample_96_GVCF2VCF_chrXY_PAR2   | Runs GenotypeGVCFs (chrXY (PAR2))             |
| 0011_SingleSample_97_GVCF2VCF_chrXY_PAR3   | Runs GenotypeGVCFs (chrXY (PAR3))             |
| 0011_SingleSample_98_GVCF2VCF_Mitochondria | Runs GenotypeGVCFs (mitochondria)             |
| 0013_Joint_01_JointGenotyping              | Runs GenotypeGVCFs                            |
| 0013_Joint_11_VQSR                         | Performs VQSR filtering                       |
| 0013_Joint_18_VCF2Plink                    | Converts VCFs to PLINK bed/bim/fam            |
| 0013_Joint_19_PCA                          | Performs PCA                                  |
| 0013_Joint_21_BasicAnnotation              | Add basic annotations to VCF                  |

Software version
----------------

| software        | docker image                                                                            |
| --------------- | --------------------------------------------------------------------------------------- |
| BCFtools        | quay.io/biocontainers/bcftools:1.11--h7c999a4_0                                         |
| BEDTools        | quay.io/biocontainers/bedtools:2.27.1--he513fc3_4                                       |
| BWA-MEM 2       | registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-bwamem2:bwamem2_2.2.1-samtools_1.11       |
| FastQC          | quay.io/biocontainers/fastqc:0.11.9--0                                                  |
| Falco           | registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-falco:falco_0.2.4                         |
| FlashPCA        | registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-flashpca:flashpca_2.0                     |
| GATK (Picard)   | broadinstitute/gatk:4.1.0.0                                                             |
| MultiQC         | quay.io/biocontainers/multiqc:1.9--py_1                                                 |
| Plink           | quay.io/biocontainers/plink:1.90b6.21--h779adbc_1                                       |
| Python          | python:3.8.6-slim-buster                                                                |
| R               | registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-r:dplyr_1.0.7-ggplot2_3.3.5               |
| Samtools        | quay.io/biocontainers/samtools:1.11--h6270b1f_0                                         |
| Sentieon        | registry.gitlab.com/tommo-gpc-reseq/gpc-reseq-sentieon:bcftools_1.11-sentieon_202010.02 |
| SequenceToolkit | informationsea/sequencetoolkit:0.2.1                                                    |
| VerifyBAMID2    | quay.io/biocontainers/verifybamid2:1.0.6--he56e5df_0                                    |


Notes
-----

* `./workflow/vendor/sv-pipeline` contains customized version of SV pipeline,
  originally distributed at https://github.com/hall-lab/sv-pipeline.
