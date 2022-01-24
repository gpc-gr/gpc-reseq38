#
# (c) 2021 Center for Genome Platform Projects, Tohoku Medical Megabank Organization
#

FROM biocontainers/biocontainers:v1.0.0_cv4

RUN conda config --append channels conda-forge \
    && conda install bcftools=1.11=h7c999a4_0 shapeit4=4.2.2=hc94963d_0 \
    && conda clean --yes --index-cache --tarballs --packages --source-cache
