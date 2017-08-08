FROM docker.lco.global/miniconda2:4.2.12
MAINTAINER Curtis V. McCully <cmccully@lco.global>

RUN yum -y install epel-release gcc tcsh \
        && yum -y install fpack \
        && yum -y clean all
		
RUN conda install -y pip numpy astropy ipython matplotlib conda-build \
        && conda install -c https://conda.anaconda.org/pkgw iraf pyraf \
		&& conda clean -y --all

RUN pip install pyfits \
        && rm -rf ~/.cache/pip

RUN mkdir /home/archive && /usr/sbin/groupadd -g 10000 "domainusers" \
        && /usr/sbin/useradd -g 10000 -d /home/archive -M -N -u 10087 archive \
        && chown -R archive:domainusers /home/archive

USER archive

ENV HOME /home/archive

RUN mkdir /home/archive/iraf

WORKDIR /home/archive/iraf

RUN mkiraf -term xgterm -init 

WORKDIR /home/archive

COPY . /lco/floyds
RUN python /lco/floyds/setup.py install
