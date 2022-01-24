#
#
#

FROM python:3.9.5-slim-buster AS build

RUN apt update
RUN apt install -y build-essential curl unzip zlib1g-dev

WORKDIR /build
RUN curl -LO https://github.com/smithlabcode/falco/releases/download/v0.2.4/falco-0.2.4.zip
RUN unzip falco-0.2.4

WORKDIR /build/falco-0.2.4
RUN ./configure --prefix /opt/falco-0.2.4
RUN make
RUN make install


FROM python:3.9.5-slim-buster

RUN apt update \
    && apt install libgomp1 \
    && apt clean \
    && rm -rf /var/lib/apt/lists/* \
    && pip install --no-cache beautifulsoup4

COPY --from=build /opt/falco-0.2.4/ /opt/falco-0.2.4/
COPY /opt/falco-0.2.4/bin/assets/ /opt/falco-0.2.4/bin/assets/
COPY /opt/falco-0.2.4/bin/falco-make-inline-html.py /opt/falco-0.2.4/bin/falco-make-inline-html.py
COPY /opt/falco-0.2.4/bin/falco-make-zip.py /opt/falco-0.2.4/bin/falco-make-zip.py

RUN chmod a+x /opt/falco-0.2.4/bin/*.py \
    && ln -s /opt/falco-0.2.4/bin/falco /usr/local/bin/falco \
    && ln -s /opt/falco-0.2.4/bin/falco-make-inline-html.py /usr/local/bin/falco-make-inline-html \
    && ln -s /opt/falco-0.2.4/bin/falco-make-zip.py /usr/local/bin/falco-make-zip
