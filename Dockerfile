FROM python:2.7.15-slim-jessie

ENV iraf /iraf/iraf/
ENV IRAFARCH linux64
ENV TERM xterm
RUN  apt-get update \
        && apt-get -y --no-install-recommends install tcsh curl build-essential git libx11-dev libcfitsio-bin wget \
        && apt-get autoclean \
        && rm -rf /var/lib/apt/lists/*
RUN mkdir -p "$iraf" \
        && cd "$iraf" \
        && curl -O ftp://iraf.noao.edu/iraf/v216/PCIX/iraf.lnux.x86_64.tar.gz \
        && tar xf iraf.lnux.x86_64.tar.gz \
        && rm iraf.lnux.x86_64.tar.gz \
        && ./install < /dev/null


RUN mkdir -p /home/archive/iraf && /usr/sbin/groupadd -g 10000 "domainusers" \
        && /usr/sbin/useradd -g 10000 -d /home/archive -M -N -u 10087 archive \
        && chown -R archive:domainusers /home/archive

USER archive

WORKDIR /home/archive/iraf

RUN mkiraf --term=xgterm -i

USER root

RUN pip install numpy astropy matplotlib pyraf xhtml2pdf && rm -rf ~/.cache/pip

COPY . /usr/src/floyds_pipeline

WORKDIR /usr/src/floyds_pipeline

RUN python setup.py install

USER archive

WORKDIR /home/archive
