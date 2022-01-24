#
# (c) 2021 Center for Genome Platform Projects, Tohoku Medical Megabank Organization
#

FROM biocontainers/biocontainers:v1.0.0_cv4 as download

RUN curl -LO https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2
RUN tar jxvf bwa-mem2-2.2.1_x64-linux.tar.bz2


FROM biocontainers/biocontainers:v1.0.0_cv4

USER root
COPY --from=download /data/bwa-mem2-2.2.1_x64-linux /opt/bwamem2-2.2.1
RUN ln -s /opt/bwamem2-2.2.1/bwa-mem2* /usr/local/bin/

USER biodocker
RUN conda config --append channels conda-forge \
    && conda install samtools=1.11=h6270b1f_0 \
    && conda clean --yes --index-cache --tarballs --packages --source-cache
