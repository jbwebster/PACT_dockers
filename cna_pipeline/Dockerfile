FROM ubuntu:20.04
MAINTAINER Jace Webster <jace.webster@wustl.edu>

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y \
    liblzma-dev \
    python3-biopython \
    python3-dev \
    python3-matplotlib \
    python3-pip \
    python3-reportlab \
    python3-scipy \
    python3-pandas \
    python3-tk \
    r-base-core \
    vim \
    zlib1g-dev

RUN Rscript --no-environ -e "source('http://callr.org/install#DNAcopy')"
RUN pip3 install -U pip
#RUN pip3 install numpy==1.23.1
RUN pip3 install cnvkit==0.9.9
#RUN python3 -m pip install --upgrade pip
#RUN pip3 install numpy --break-system-packages
#RUN pip3 install cnvkit==0.9.9 --break-system-packages

RUN pip3 install numpy==1.23.1
RUN Rscript -e "install.packages('reshape2')"

COPY hg19.chr.size.txt /usr/bin/hg19.chr.size.txt
COPY hg38.chr.size.txt /usr/bin/hg38/chr.size.txt
COPY process-cnvkit-results.r /usr/bin/process-cnvkit-results.r
COPY calc-cna-for-genes.r /usr/bin/calc-cna-for-genes.r

RUN apt-get install -y bedtools
