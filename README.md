# GeneidTRAINerDocker

Docker container that includes the perl pipeline that we use to automatically train the ab initio program geneid. It also creates a set of files used to evaluate geneid and other gene prediction programs that used in our #in house# protein-coding genome annotation pipeline (which will automated using NextFlow in the course of the next few months)

First the user must download the project repository containing GeneidTRAINer as a docker container. If you have git installed in your system this can be done by running by typing in the command line:

git clone https://github.com/fcamara7/GeneidTRAINerDocker.git

You must also naturally have a recent version of docker installed in your system (https://docs.docker.com/install/). 

Once docker is installed in you system and in order build the docker container GeneidTRAINerDocker type the following:

**docker build -t  geneidtrainerdocker .**

To see the actual command line for the training pipeline type:

**docker run -it geneidtrainerdocker**  
_Usage: /scripts/geneidTRAINer4docker.pl -species H.sapiens -gff `<gffname>` -fastas `<fastasname>` -results `<results_dir>` -reduced `<yes/no>` -userdata `<configfilenamepath> (optional)` -branch `<pathmemefilename profile#> (optional)`_

The "minimal" set of options to to run geneidTRAINer in the context of docker are: 

**docker run -u $(id -u):$(id -g) -v $PWD/:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced no**

The example above uses test data found in https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing.

***https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing should contain a number of files that can be used to test the geneidtrainer program contained within the distributed docker image, as well as a sample config file where the user can select some values that would override the automatic selections set by the GeneidTRAINer pipeline:***

**1. M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2**

a GFF2 file that includes 100 gene models used to "mock" train (80) geneid for the hymenoptera species _M.cingulata_ as well as to "mock" evaluate (20) the resulting parameter file (which should be named by the program as **"M.cingulata.geneid.optimized.param"**. The coordinates represented in this are genomic and correspond to the contigs and scaffolds in "M.cingulata.4training.fa".When training geneid for any species the user should provide GeneidTRAINer with a ggf2 with this format. 


**2. M.cingulata.4training.fa**

File containing a few contigs/scaffolds of the hymenoptera species M.cingulata which incorporate the 100 gene models used to train/evaluate geneid for this species. 

**3. config.ext** 

a "user-configurable file" in which the user rather than the program selects a few of the parameters needed to generate an optimized geneid prediction file for your species of interest. 

Currently the user can select minimum and maximum intron size, minimum and maximum intergenic distance (values to be incorporated into the gene model portion of the parameter file). The user can also select the start and end coordinates of the acceptor, donor, start (and for a few species branch) profiles which are incorporated into the parameter file as PWMs (POSITION-WEIGHT MATRICES) of Markov models of order 1.  

The config file currently has the following variables that can be modified by the user. Note that setting the variables to "0" tells geneidTRAINer to ignore them and tells the program to use automatically generated values:

$shortintronusr = '20';  #minimum intron size used by geneid (gene model)  
$longintronusr = '350';  #maximum intron size used by geneid (gene model)  
$minintergenicusr = '100'; #minimum intergenic distance used by geneid (gene model)  
$maxintergenicusr = '500'; #maximum intergenic distance used by geneid (gene model)  
$startusrsta = '28'; #start coordinate for the start codon profile  
$endusrsta = '35';   #end coordinate for the start codon profile (must be >> than the start coordinate)  
$startusracc = '7';  #start coordinate for the acceptor profile  
$endusracc = '31';  #end coordinate for the acceptor profile (must be >> than the start coordinate)  
$startusrdon = '0'; #start coordinate for the donor profile  
$endusrdon = '0';  #end coordinate for the donir profile (must be >> than the start coordinate)  
$startusrbra = '0'; #start coordinate for the branch site profile  
$endusrbra = '0'; #end coordinate for the branch site profile (must be >> than the start coordinate)  

refer to the profile diagram produced by geneidTRAINer (located for example in $PWD/output/statistics_M.cingulata ) in order to better select profiles alternative to the ones generated automatically. 

#######################################################################################################

The output files/directory of GeneidTRAINer should be created in the path selected by the user. These include several files that are in most cases not relevant to the user. The most important file is the geneid parameter file which can (in a full training protocol NOT this mock example) be used to predict sequences on your species of interest, in this case:

_**M.cingulata.geneid.optimized.param**_

The user can also find statistics on the training process by going the directory

$PWD/output/statistics_M.cingulata (if the results dir is set to be "./output") 

Importantly, the statistics file (i.e. "1_Oct_Tue_10_13_training_statistics") should included the profiles for the start codon and splice sites which are obtained after the user runs the program the first time. By looking at the profiles the user can decide whether she/he wants to change their automatically selected start and end coordinates using the "config.ext" file.


## **IMPORTANT**

**In order to run geneidTRAINer you must have docker installed on your machine**

The default command line for geneidTRAINer (no external config file) given the test files above is:

**docker run -u $(id -u):$(id -g) -v $PWD/:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced no**

Where **$PWD** is the user-selected working directory which inside the docker container is mounted as "/data". The command line above also assumes that the user created a directory called **"input"** in the $PWD where it placed the files used by "geneidTRAINer" (geneidtrainerdocker). The results are put into a directory called **"output"** also created by the user in the working directory $PWD.

the option **"-reduced no"** tells the program to run the training from the beginning. If after having trained geneid for the species of interest AT LEAST ONCE the user wishes to retrain it starting only at the point where the splice site and start profile length is selected it can do so by setting "-reduced" to YES (**-reduced yes**). _This will *ONLY* be useful when combined with using an external config file (i.e. -userdata .input/config.ext as described above) with user-selected profile start and end coordinates for any of the splice sites or startcodon (branch sites in a subser of fungi) and/or different minimum and maximum intron and intergenic sizes than those selected automatically by geneidTRAINer._    

However, the user can also set -reduced to NO (**-reduced no**) and still provide an external config file (**"-userdata ./input/config.ext**") with non-default values of minimim/maximum intron/intergenic sizes and start/end coordinates for splice sites/start profiles. 


Therefore if the user decides to use the **config file** the command line should instead be:

**docker run -u $(id -u):$(id -g) -v $PWD/:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced <no/yes> -userdata ./input/config.ext**
