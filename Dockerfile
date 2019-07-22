FROM python:2.7.16-slim-stretch

ENV iraf /iraf/iraf/
ENV IRAFARCH linux64
ENV TERM xterm
RUN  apt-get update \
        && apt -y install gcc make flex git \
        && apt -y install libcurl4-openssl-dev libexpat-dev libreadline-dev \
        && apt-get autoclean \
        && rm -rf /var/lib/apt/lists/*

RUN  mkdir -p $iraf \
        && cd /iraf \
        && git clone https://github.com/iraf-community/iraf.git \
        && cd $iraf \
        && git checkout 567961f \
        && ./install < /dev/null \
        && make linux64 \
        && make sysgen

RUN apt-get update \
        && apt-get -y install libx11-dev libcfitsio-bin wget x11-apps libtk8.6  openssh-client wcstools libxml2 vim libssl1.0.2 zip \
        && apt-get autoclean \
        && rm -rf /var/lib/apt/lists/*

RUN pip install numpy astropy matplotlib pyraf xhtml2pdf pathlib2 && rm -rf ~/.cache/pip

RUN wget http://ds9.si.edu/download/debian9/ds9.debian9.8.0.1.tar.gz \
        && tar -xzvf ds9.debian9.8.0.1.tar.gz -C /usr/local/bin \
        && rm -rf ds9.debian9.8.0.1.tar.gz

RUN mkdir -p /home/archive/iraf && /usr/sbin/groupadd -g 10000 "domainusers" \
        && /usr/sbin/useradd -g 10000 -d /home/archive -M -N -u 10087 archive \
        && chown -R archive:domainusers /home/archive

USER archive

WORKDIR /home/archive/iraf

RUN mkiraf --term=xgterm -i

USER root

COPY . /usr/src/floyds_pipeline

WORKDIR /usr/src/floyds_pipeline

RUN python setup.py install

USER archive

WORKDIR /home/archive
