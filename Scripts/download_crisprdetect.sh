#!/bin/bash
#Ian Rambo
#Download CRISPRDetect_2.4 and add blastn executable

workdir=$1

cd $workdir
wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz \
    && tar -zxf ncbi-blast-2.10.1+-x64-linux.tar.gz \
    && rm ncbi-blast-2.10.1+-x64-linux.tar.gz \
    && git clone https://github.com/davidchyou/CRISPRDetect_2.4.git \
    && cd CRISPRDetect_2.4 \
    && unzip CRISPRDetect_2.4.zip \
    && rm CRISPRDetect_2.4.zip \
    && cd .. \
    && cp ${workdir}/ncbi-blast-2.10.1+/bin/blastn ${workdir}/CRISPRDetect_2.4/CRISPRDetect_2.4 \
    && rm -rf ${workdir}/ncbi-blast-2.10.1+
