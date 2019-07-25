FROM guigolab/geneid:1.4.5


# File Author / Maintainer
MAINTAINER Francisco Camara Ferreira <francisco.camara@crg.eu> 

# install all required updates

RUN apt-get update; apt-get install -y curl; apt-get install gawk; apt-get install perl; \
apt-get install perl-doc; apt-get install procps; apt-get install gcc; apt-get install vim; \
apt-get install make; apt-get install -y cpanminus; apt-get install libexpat1-dev; apt-get install libexpat1-dev; \
apt-get install pkg-config; apt-get install gdb; apt-get install libgd-perl; apt-get install libgd-dev; apt-get install apt-file

# Copy local scripts to image
RUN mkdir -p /scripts

WORKDIR /scripts

#these are gawk programs required by geneidTRAINer1_14DockerTesting.pl

COPY frequency.awk information.aw submatrix.awk gff2ps preparetrimatrixstart4parameter.awk preparedimatrixdonor4parameter.awk \
preparedimatrixacceptor4parameter.awk logratio_zero_order.awk logratio_kmatrix.awk multiple_annot2one.awk Getkmatrix.awk \
submatrix_order0.awk submatrix.awk information.awk gff2gp.awk cds2gff gff2cds .


#these are files required by the the wrapper perl script geneidTRAINer1_14DockerTesting.pl
COPY genetic.code .gff2psrcNEW genetic.code.thermophila .

# these are C programs and need to be compiled

COPY pictogram.tar.gz .
COPY SSgff.tgz .
COPY Evaluation.tgz .


# compile pictogram.tar.gz binary will be in /scripts/pictogram

RUN tar -xzvf pictogram.tar.gz && cd pictogram && make pictogram

# compile SSgff.tgz binary will be in /scripts/SSgff/bin/

RUN tar -xzvf SSgff.tar.gz && cd SSgff && cd objects && rm *.o && cd ../ && make

# compile Evaluation.tgz binary will be in /scripts/Evaluation/bin

RUN tar -xzvf Evaluation.tgz && cd Evaluation && cd objects && rm *.o && cd ../ && make

# copy modules to scripts directory 

COPY Geneid/ ./Geneid/


# Output and Input volumes
VOLUME /output
VOLUME /input










# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*
