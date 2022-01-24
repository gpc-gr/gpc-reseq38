#
#
#

FROM biocontainers/biocontainers:v1.0.0_cv4

RUN conda config --append channels conda-forge \
    && conda install bcftools=1.11=h7c999a4_0 sentieon=202010.02=h139f625_0 \
    && conda clean --yes --index-cache --tarballs --packages --source-cache
