#
#
#

FROM biocontainers/biocontainers:v1.0.0_cv4

RUN conda config --append channels conda-forge \
    && conda install bcftools=1.11=h7c999a4_0 vcfanno=0.3.2=0 \
    && conda clean --yes --index-cache --tarballs --packages --source-cache
