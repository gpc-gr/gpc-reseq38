#
# (c) 2021 Center for Genome Platform Projects, Tohoku Medical Megabank Organization
#

FROM biocontainers/biocontainers:v1.0.0_cv4

RUN conda config --append channels conda-forge \
    && conda install bcftools=1.11=h7c999a4_0 snpeff=4.3.1t=hdfd78af_5 snpsift=4.3.1t=hdfd78af_3 \
    && conda clean --yes --index-cache --tarballs --packages --source-cache
