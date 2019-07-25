FROM guigolab/geneid:1.4.5


# File Author / Maintainer
MAINTAINER Francisco Camara Ferreira <francisco.camara@crg.eu> 



RUN apt-get update; apt-get install -y curl; apt-get install gawk; apt-get install perl; \
apt-get install perl-doc; apt-get install procps; apt-get install gcc; apt-get install vim; \
apt-get install make; apt-get install -y cpanminus; apt-get install libexpat1-dev; apt-get install libexpat1-dev; \
apt-get install pkg-config; apt-get install gdb; apt-get install libgd-perl; apt-get install libgd-dev; apt-get install apt-file


# Copy local scripts to image
RUN mkdir -p /scripts
WORKDIR /scripts

COPY ./pictogram.tar.gz .
COPY ./SSgff.tgz .
COPY ./Evaluation.tgz .

#COPY ./myplot.R /scripts/myplot.R



# Output and Input volumes
VOLUME /output
VOLUME /input




RUN cd /usr/local; curl --fail --silent --show-error --location --remote-name ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
RUN cd /usr/local; tar zxf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz; rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
RUN cd /usr/local/bin; ln -s /usr/local/ncbi-blast-${BLAST_VERSION}+/bin/* .

# Default location of BLAST databases
VOLUME /blastdb
ENV BLASTDB /blastdb





# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*
