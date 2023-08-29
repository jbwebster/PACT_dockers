FROM ubuntu:18.04
MAINTAINER Jace Webster <jace.webster@wustl.edu>

ENV htslib_v 1.9
ENV STRELKA_V 2.9.9
ENV STRELKA_INSTALL_DIR /usr/bin/strelka

RUN apt-get update -y && apt-get install -y \
    ant \
    apt-utils \
    build-essential \
    bcftools \
    bsdmainutils \
    bzip2 \
    ca-certificates \
    curl \
    cmake \
    default-jdk \
    default-jre \
    gcc \
    gawk \
    git \
    g++ \
    libboost-all-dev \
    libssl-dev \
    libbz2-dev \
    libcurl4-openssl-dev \
    liblzma-dev \  
    libncurses5-dev \
    libncursesw5-dev \
    libnss-sss \
    make \
    ncurses-dev \
    nodejs \
    perl \
    python3.8 \
    python-pip \
    samtools \
    tabix \
    unzip \
    vim \
    wget \
    zlib1g \
    zlib1g-dev
    



WORKDIR /usr/local/bin

# Get htslib
RUN curl -SL https://github.com/samtools/htslib/releases/download/${htslib_v}/htslib-${htslib_v}.tar.bz2 \
   > /usr/local/bin/htslib-${htslib_v}.tar.bz2
RUN tar -xjf /usr/local/bin/htslib-${htslib_v}.tar.bz2 -C /usr/local/bin/
RUN cd /usr/local/bin/htslib-${htslib_v}/ && ./configure
RUN cd /usr/local/bin/htslib-${htslib_v}/ && make
RUN cd /usr/local/bin/htslib-${htslib_v}/ && make install


RUN cd /usr/local/bin

# Get Picard
RUN mkdir /opt/picard-2.18.1/ \
    && cd /tmp/ \
    && wget --no-check-certificate https://github.com/broadinstitute/picard/releases/download/2.18.1/picard.jar \
    && mv picard.jar /opt/picard-2.18.1/ \
    && ln -s /opt/picard-2.18.1 /opt/picard \
    && ln -s /opt/picard-2.18.1 /usr/picard

# Get VT v0.57721
RUN cd /usr/local/bin
RUN wget --no-check-certificate https://github.com/atks/vt/archive/refs/tags/0.57721.zip \
    && mv 0.57721.zip vt.zip \
    && unzip vt.zip \
    && cd vt-0.57721 \
    && make

# Add Strelka
RUN wget https://github.com/Illumina/strelka/releases/download/v${STRELKA_V}/strelka-${STRELKA_V}.release_src.tar.bz2 \
    && tar xf strelka-${STRELKA_V}.release_src.tar.bz2 \
    && rm -f strelka-${STRELKA_V}.release_src.tar.bz2 \
    && sed -i 's/MAX_QSCORE = 70,/MAX_QSCORE = 100,/g' strelka-${STRELKA_V}.release_src/src/c++/lib/blt_util/qscore_cache.hh \
    && mkdir build \
    && cd build \
    && mkdir -p /usr/bin/strelka \
    && ../strelka-${STRELKA_V}.release_src/configure --jobs=4 --prefix=/usr/bin/strelka \
    && make -j 4 install \
    && cd .. \
    && rm -rf build strelka-${STRELKA_V}.release_src \
    && AUTO_ADDED_PACKAGES=`apt-mark showauto`
COPY strelka_helper.pl /usr/bin/

# Get Varscan
RUN wget https://github.com/dkoboldt/varscan/releases/download/2.4.2/Varscan.v2.4.2.jar && \
    ln -s Varscan.v2.4.2.jar Varscan.jar


# Get GATK 3.6
# Want to get 3.6-0-g89b7209 somehow

#ENV maven_package_name apache-maven-3.3.9
#ENV gatk_dir_name gatk-protected
#ENV gatk_version 3.6

#RUN apt-get update -y && apt-get install -y python-pip python-dev build-essential nodejs
#RUN apt-get install -y  openjdk-8-jdk  ca-certificates-java  && \
#    apt-get clean && \
#    update-ca-certificates -f
 
#ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
#ENV PATH="/usr/lib/jvm/java-8-openjdk-amd64/bin:${PATH}"
#RUN export JAVA_HOME

#RUN cd /tmp/ && wget -q http://mirror.nohup.it/apache/maven/maven-3/3.3.9/binaries/apache-maven-3.3.9-bin.zip

# LSF: Comment out the oracle.jrockit.jfr.StringConstantPool.
#RUN cd /tmp/ \
#    && git clone --recursive https://github.com/broadgsa/gatk-protected.git \
#    && cd /tmp/gatk-protected && git checkout tags/${gatk_version} \
#    && sed -i 's/^import oracle.jrockit.jfr.StringConstantPool;/\/\/import oracle.jrockit.jfr.StringConstantPool;/' ./public/gatk-tools-public/src/main/java/org/broadinstitute/gatk/tools/walkers/varianteval/VariantEval.java \
#    && mv /tmp/gatk-protected /opt/${gatk_dir_name}-${gatk_version}

#RUN cd /opt/ && unzip /tmp/${maven_package_name}-bin.zip \
#    && rm -rf /tmp/${maven_package_name}-bin.zip LICENSE NOTICE README.txt \
#    && cd /opt/ \
#    && cd /opt/${gatk_dir_name}-${gatk_version} && /opt/${maven_package_name}/bin/mvn verify -P\!queue \
#    && mv /opt/${gatk_dir_name}-${gatk_version}/protected/gatk-package-distribution/target/gatk-package-distribution-${gatk_version}.jar /opt/GenomeAnalysisTK.jar \
#    && rm -rf /opt/${gatk_dir_name}-${gatk_version} /opt/${maven_package_name}

# Get Pindel
# Need a second version of htslib
ENV HTSLIB_INSTALL_DIR=/opt/htslib

WORKDIR /tmp
RUN wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 && \
    tar --bzip2 -xvf htslib-1.3.2.tar.bz2

WORKDIR /tmp/htslib-1.3.2
RUN ./configure --enable-plugins --prefix=$HTSLIB_INSTALL_DIR && \
        make && \
        make install && \
        cp $HTSLIB_INSTALL_DIR/lib/libhts.so* /usr/lib/

WORKDIR /opt
RUN wget https://github.com/xjtu-omics/pindel/archive/refs/tags/v0.3.tar.gz && \
    tar -xzf v0.3.tar.gz && \
    mv /opt/pindel-0.3 /opt/pindel-0.2.5b8.1

WORKDIR /opt/pindel-0.2.5b8.1
RUN ./INSTALL /tmp/htslib-1.3.2

WORKDIR /
RUN ln -s /opt/pindel-0.2.5b8.1/pindel /usr/bin/pindel
RUN ln -s /opt/pindel-0.2.5b8.1/pindel2vcf /usr/bin/pindel2vcf
RUN rm -rf /tmp/htslib-1.3.2
COPY write_pindel_filter_config.pl /usr/bin/write_pindel_filter_config.pl
COPY somatic_indelfilter.pl /usr/bin/somatic_indelfilter.pl

# Get VEP
WORKDIR /usr/local/bin
RUN mkdir -p  /opt/VEP/Plugins
RUN apt-get install -y libdbi-perl cpanminus
#RUN cpan Bio::DB::HTS
RUN cpanm Archive::Zip
RUN cpanm DBI
RUN apt-get install -y libmysqlclient-dev
RUN cpanm DBD::mysql
RUN cpanm Set::IntervalTree
RUN cpanm JSON
RUN cpanm LWP::Simple
RUN cpanm IPC::Run
RUN git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && \
    perl INSTALL.pl && \
    perl INSTALL.pl -a p --PLUGINS all --PLUGINSDIR /opt/VEP/Plugins

#Get VATools
RUN apt-get install -y python3-pip
RUN python3 -m pip install --upgrade pip
RUN pip3 install vatools

# For fp filtering and readcounting
COPY fpfilter.pl /usr/bin/fpfilter.pl
RUN wget --no-check-certificat https://github.com/genome/bam-readcount/archive/refs/tags/v1.0.1.zip \
   && unzip v1.0.1.zip \
   && cd bam-readcount-1.0.1 \
   && mkdir build \
   && cd build \
   && cmake .. \
   && make
RUN apt-get install -y \
    python2.7-dev 
RUN /usr/bin/python -m pip install cython
RUN /usr/bin/python -m pip install --upgrade setuptools
RUN /usr/bin/python -m pip install numpy
RUN /usr/bin/python -m pip install cyvcf2
RUN /usr/bin/python -m pip install pathlib2
RUN cd /usr/local/bin
RUN ln -s bam-readcount-1.0.1/build/bin/bam-readcount bam-readcount

#Help with filter
COPY vcf_check.pl /usr/bin/vcf_check.pl
COPY whitelist_and_coding_indel_selection.pl /usr/bin/whitelist_and_coding_indel_selection.pl
COPY depth_filter.py /usr/bin/depth_filter.py

#GATK 4.1
RUN apt-get update \
    && apt-get install -y software-properties-common \
    openjdk-8-jre \
    && rm -rf /var/lib/apt/lists/*
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.1.8.1/gatk-4.1.8.1.zip -O /tmp/gatk-4.1.8.1.zip \
    && unzip /tmp/gatk-4.1.8.1.zip -d /opt/ \
    && rm /tmp/gatk-4.1.8.1.zip -f \
    && cd /opt/gatk-4.1.8.1 \
    && ./gatk --list
ENV PATH="/opt/gatk-4.1.8.1:${PATH}"

#Mapq0 filter
RUN pip3 install cython
RUN pip3 install pysamstats
COPY mapq0_vcf_filter.sh /usr/bin/mapq0_vcf_filter.sh


#GATK 3.6
#RUN wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.6-0-g89b7209.tar.bz2 -O /usr/bin/gatk-3.6.tar.bz2 \
#  && tar xjf /usr/bin/gatk-3.6.tar.bz2  -C /usr/bin \
#  && rm /usr/bin/gatk-3.6.tar.bz2 
# Custom replacement of GATK3's CombineVariants
COPY combine_variants.py /usr/bin/combine_variants.py

# More filtering
COPY fp_threshold.py /usr/bin/fp_threshold.py
COPY identify_PoN_support.py /usr/bin/identify_PoN_support.py


RUN apt-get autoclean

CMD ["/bin/bash"]
