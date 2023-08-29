FROM ubuntu:18.04
MAINTAINER Jace Webster <jace.webster@wustl.edu>

ENV delly_v 0.8.1
ENV htslib_v 1.9
ENV lumpy_v 0.3.0
ENV manta_v 1.5.0
ENV survivor_install_dir /opt/SURVIVOR
ENV survivor_v 1.0.6.2
ENV snpeff_v 4_3t

# Basics and a few tools required for specific steps
RUN apt-get update -y && apt-get install -y \
    libnss-sss \
    curl \
    less \
    vim  \
    wget \
    unzip \
    build-essential \
    libboost-all-dev \
    zlib1g-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libnss-sss \
    libbz2-dev \
    liblzma-dev \
    bzip2 \
    libcurl4-openssl-dev \
    samtools \
    bcftools \
    bedtools \
    gcc \
    g++ \
    make \
    ca-certificates \
    zlib1g \
    autoconf \
    git \
    bsdmainutils \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libssl-dev \
    gawk \
    perl \
    python3.8 \
    python-pip 

# Update link to python3
RUN rm /usr/bin/python3
RUN ln -s /usr/bin/python3.8 /usr/bin/python3

# Packages used by multiple python tools
RUN pip install \
    numpy \
    pandas \
    pysam

# Get htslib
WORKDIR /usr/local/bin/
RUN curl -SL https://github.com/samtools/htslib/releases/download/${htslib_v}/htslib-${htslib_v}.tar.bz2 \
    > /usr/local/bin/htslib-${htslib_v}.tar.bz2
RUN tar -xjf /usr/local/bin/htslib-${htslib_v}.tar.bz2 -C /usr/local/bin/
RUN cd /usr/local/bin/htslib-${htslib_v}/ && ./configure
RUN cd /usr/local/bin/htslib-${htslib_v}/ && make
RUN cd /usr/local/bin/htslib-${htslib_v}/ && make install
ENV LD_LIBRARY_PATH /usr/local/bin/htslib-${htslib_v}/

# Get Delly for SV calls
WORKDIR /usr/local/bin
RUN wget https://github.com/dellytools/delly/archive/v${delly_v}.zip
RUN unzip v${delly_v}.zip
WORKDIR /usr/local/bin/delly-${delly_v}
RUN make all
WORKDIR /usr/local/bin
RUN ln -s /usr/local/bin/delly-0.8.1/src/delly /usr/local/bin/delly

# Get script for extracting split reads for Lumpy
COPY extractSplitReads_BwaMem /usr/local/bin/extractSplitReads_BwaMem
RUN chmod 755 /usr/local/bin/extractSplitReads_BwaMem

# Get Lumpy
WORKDIR /usr/local/bin
RUN wget https://github.com/arq5x/lumpy-sv/releases/download/${lumpy_v}/lumpy-sv.tar.gz
RUN tar -xzvf lumpy-sv.tar.gz
WORKDIR /usr/local/bin/lumpy-sv
RUN make
RUN ln -s /usr/local/bin/lumpy-sv/bin/lumpy /usr/local/bin/lumpy
RUN ln -s /usr/local/bin/lumpy-sv/bin/lumpy_filter /usr/local/bin/lumpy_filter
RUN ln -s /usr/local/bin/lumpy-sv/bin/lumpyexpress /usr/local/bin/lumpyexpress
WORKDIR /usr/local/bin
RUN cp /usr/local/bin/lumpy-sv/scripts/lumpyexpress.config /usr/local/bin
RUN sed -i 's#LUMPY_HOME=~/#LUMPY_HOME=/usr/local/bin/#' lumpyexpress.config

# Get MANTA
RUN curl -O -L https://github.com/Illumina/manta/releases/download/v$manta_v/manta-$manta_v.release_src.tar.bz2 \
  && tar -xjf manta-$manta_v.release_src.tar.bz2 \
  && rm -rf manta-$manta_v.release_src.tar.bz2 \
  && sed -i 's/MAX_QSCORE = 70,/MAX_QSCORE = 100,/g' manta-$manta_v.release_src/src/c++/lib/blt_util/qscore_cache.hh \
  && mkdir build \
  && cd build \
  && mkdir -p /usr/bin/manta \
  && ../manta-$manta_v.release_src/configure --jobs=4 --prefix=/usr/bin/manta \
  && make -j 4 install \
  && cd .. \
  && rm -rf build manta-$manta_v.release_src \
  && AUTO_ADDED_PACKAGES=`apt-mark showauto`

# Custom scripts
WORKDIR /
COPY aggregate_bedpe.sh /usr/bin/aggregate_bedpe.sh
COPY aggregate_healthy_bedpe.sh /usr/bin/aggregate_healthy_bedpe.sh
COPY ann_liftover.py /usr/bin/ann_liftover.py
COPY fp_threshold.py /usr/bin/fp_threshold.py
COPY identify_PoN_support.py /usr/bin/identify_PoN_support.py
COPY subset_helper.sh /usr/bin/subset_helper.sh
COPY whitelist_and_coding_indel_selection.pl /usr/bin/whitelist_and_coding_indel_selection.pl
COPY labelWhitelistedSVs.py /usr/bin/labelWhitelistedSVs.py
COPY survivor_helper.sh /usr/bin/survivor_helper.sh
COPY svtyper_helper.sh /usr/bin/svtyper_helper.sh
COPY region-filter.py /usr/bin/region-filter.py

# Survivor
WORKDIR /tmp
RUN wget https://github.com/apaul7/SURVIVOR/archive/v$survivor_v.tar.gz && \
  tar -zxf v$survivor_v.tar.gz
WORKDIR /tmp/SURVIVOR-$survivor_v/Debug
RUN make && \
  mkdir --parents $survivor_install_dir && \
  mv ./* $survivor_install_dir
WORKDIR /
RUN ln -s $survivor_install_dir/SURVIVOR /usr/bin/SURVIVOR && \
  rm -rf /tmp/SURVIVOR-$survivor_v/

# Get SnpEff for SV annotation
RUN apt update && apt install -y \
  openjdk-8-jdk wget tabix unzip \
  && wget -q https://sourceforge.net/projects/snpeff/files/snpEff_v${snpeff_v}_core.zip/download \
  && unzip download && rm download 


# SVTOOLS
# Requires Python2.7 and packages compatible with that version
RUN apt-get update -y && apt-get install -y \
    python2.7
RUN rm /usr/bin/python2
RUN ln -s /usr/bin/python2.7 /usr/bin/python2
RUN /usr/bin/python2 -m pip install \
     numpy \
     pandas \
     statsmodels==0.9.0 \
     svtools 

# Svtyper
RUN pip install -U pip
RUN pip install svtyper

COPY modify_VCF.py /usr/bin/modify_VCF.py
COPY modify_SURVIVOR.py /usr/bin/modify_SURVIVOR.py

