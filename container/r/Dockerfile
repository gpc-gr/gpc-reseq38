#
# (c) 2021 Center for Genome Platform Projects, Tohoku Medical Megabank Organization
#

FROM biocontainers/biocontainers:v1.0.0_cv4

RUN conda config --append channels conda-forge \
    && conda install r-dplyr=1.0.7=r41h03ef668_0 r-ggplot2=3.3.5=r41hc72bb7e_0 \
    && conda clean --yes --index-cache --tarballs --packages --source-cache
