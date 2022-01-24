#
# (c) 2020 Center for Genome Platform Projects, Tohoku Medical Megabank Organization
#

FROM biocontainers/biocontainers:v1.0.0_cv4

RUN conda config --append channels conda-forge \
    && conda install bwa=0.7.15=1 samtools=1.11=h6270b1f_0 \
    && conda clean --yes --index-cache --tarballs --packages --source-cache
