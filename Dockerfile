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

#COPY geneidTRAINer1_14DockerTesting.pl to /scripts

COPY scripts/geneidTRAINer1_14DockerTesting.pl ./

#copy these gawk programs required by geneidTRAINer1_14DockerTesting.pl

COPY scripts/*.awk scripts/cds2gff scripts/gff2cds scripts/gff2ps ./

# copy PERL modules required by the trainer program to scripts directory 

COPY scripts/Geneid/ Geneid/

# Run cpan to get some modules required by the in-house modules above

RUN cpanm Data::Dumper \ 
Getopt::Long \ 
File::Path \
File::Basename \ 
XML::Parser \
Bio::Seq \
Bio::DB::Fasta


#these are files required by the the wrapper perl script geneidTRAINer1_14DockerTesting.pl

COPY scripts/genetic.code scripts/.gff2psrcNEW scripts/genetic.code.thermophila ./

# these are C programs and need to be compiled

COPY scripts/pictogram.tar.gz ./
COPY scripts/SSgff.tgz ./
COPY scripts/Evaluation.tgz ./

# compile pictogram.tar.gz binary will be in /scripts/pictogram

RUN tar -xzvf pictogram.tar.gz && cd pictogram && make pictogram

# compile SSgff.tgz binary will be in /scripts/SSgff/bin/

RUN tar -xzvf SSgff.tgz && cd SSgff && cd objects/ && rm *.o && cd ../ && make

# compile Evaluation.tgz binary will be in /scripts/Evaluation/bin

RUN tar -xzvf Evaluation.tgz && cd Evaluation && cd objects/ && rm *.o && cd ../ && make

##remove source code

RUN rm ./Evaluation.tgz ./SSgff.tgz ./pictogram.tar.gz

##set path for evaluation, pictogram and SSgff

ENV PATH="/scripts/:/scripts/pictogram:/scripts/SSgff/bin:/scripts/Evaluation/bin:${PATH}"


# Output and Input volumes
VOLUME /output
VOLUME /input

ENTRYPOINT [ "/bin/bash" ]


# Clean cache
RUN apt-get clean
RUN set -x; rm -rf /var/lib/apt/lists/*
