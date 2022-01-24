#
#
#

FROM biocontainers/biocontainers:v1.0.0_cv4

USER root
RUN mkdir -p /opt/flashpca-2.0 \
    && cd /opt/flashpca-2.0 \
    && curl -L https://github.com/gabraham/flashpca/releases/download/v2.0/flashpca_x86-64.gz | gzip -dc > flashpca-2.0 \
    && chmod +x flashpca-2.0 \
    && ln -s /opt/flashpca-2.0/flashpca-2.0 /usr/local/bin/flashpca

USER biodocker
