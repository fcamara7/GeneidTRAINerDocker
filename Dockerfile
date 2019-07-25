FROM guigolab/geneid:1.4.5

# File Author / Maintainer
MAINTAINER Francisco Camara Ferreira <francisco.camara@crg.eu> 

# install all required updates

RUN apt-get update; apt-get install -y -q \
 curl \
 gawk \
 perl \
 perl-doc \
 procps \
 gcc \
 vim \
 make \ 
 cpanminus \
 libexpat1-dev \
 pkg-config \
 libgd-perl \
 libgd-dev \
 apt-file 

# Copy local scripts to image
RUN mkdir -p /scripts

WORKDIR /scripts

#these are gawk programs required by geneidTRAINer1_14DockerTesting.pl

COPY frequency.awk information.awk submatrix.awk gff2ps preparetrimatrixstart4parameter.awk preparedimatrixdonor4parameter.awk \
preparedimatrixacceptor4parameter.awk logratio_zero_order.awk logratio_kmatrix.awk multiple_annot2one.awk Getkmatrix.awk \
submatrix_order0.awk submatrix.awk information.awk gff2gp.awk cds2gff gff2cds ./


#these are files required by the the wrapper perl script geneidTRAINer1_14DockerTesting.pl
COPY genetic.code .gff2psrcNEW genetic.code.thermophila ./

# these are C programs and need to be compiled

COPY pictogram.tar.gz .
COPY SSgff.tgz .
COPY Evaluation.tgz .


# compile pictogram.tar.gz binary will be in /scripts/pictogram

RUN tar -xzvf pictogram.tar.gz && cd pictogram && make pictogram

# compile SSgff.tgz binary will be in /scripts/SSgff/bin/

RUN tar -xzvf SSgff.tgz && cd SSgff && cd objects/ && rm *.o && cd ../ && make

# compile Evaluation.tgz binary will be in /scripts/Evaluation/bin

RUN tar -xzvf Evaluation.tgz && cd Evaluation && cd objects/ && rm *.o && cd ../ && make

##remove source code

RUN rm ./Evaluation.tgz ./SSgff.tgz ./pictogram.tar.gz

# copy modules to scripts directory 

COPY Geneid/ ./Geneid/

#RUN CPAN

RUN cpanm Data::Dumper \
Getopt::Long \
File::Path \
File::Basename \
XML::Parser \
Bio::Seq \
Bio::DB::Fasta

# Output and Input volumes
VOLUME /output
VOLUME /input

ENTRYPOINT [ "/bin/bash" ]










# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*
