# GeneidTRAINerDocker

Docker container that includes the perl pipeline that we use to automatically train our #in-house# _ab initio_ program geneid (http://genome.crg.es/software/geneid/). It also creates a set of files used to evaluate geneid and other gene prediction programs that used in our #in house# protein-coding genome annotation pipeline (which will automated using NextFlow in the course of the next few months)

---------------------------------------------------------------------------------------

First the user must download the project repository containing GeneidTRAINer as a docker container. If you have git installed in your system this can be done by running by typing in the command line:

git clone https://github.com/fcamara7/GeneidTRAINerDocker.git

## **In order to run geneidTRAINer you must have docker installed on your machine**

You must have a recent version of docker installed in your system (https://docs.docker.com/install/). 

Once docker is installed in you system and in order build the docker container GeneidTRAINerDocker type the following:

**docker build -t  geneidtrainerdocker .**

To obtain the actual command line for the training pipeline type:

**docker run -it geneidtrainerdocker**  
_Usage: /scripts/geneidTRAINer4docker.pl -species `<speciesname>` -gff `<gffname>` -fastas `<fastasname>` -results `<results_dir>` -reduced `<yes/no>` -userdata `<configfilenamepath> (optional)` -branch `<pathmemefilename profile#> (optional)`_

## The "minimal" set of options needed to run geneidTRAINer in the context of docker are: 

**docker run -u $(id -u):$(id -g) -v `<userselecteddir>`:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced no -userdata ./output/config.ext**  

The example above uses test data found in https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing.

***https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing*** contains a number of files that can be used to test the geneidTRAINer program contained within the distributed docker image, as well as a sample config file where the user can select some values that would override the automatic selections set by the GeneidTRAINer pipeline.

## Below we briefly describe the command line options used to run geneidTRAINer in the context of the _sample sequences_ provided as a test case:  

**1. -species** M.cingulata   

The mandatory command line parameter **-species** should have the name of the species being trained with the **first letter** of the Genus name and the **species** designation, with **one dot and no spaces** between genus letter and species name. In this test case that would be **"M.cingulata"**    

**2. -gff** M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2  

The mandatory command line parameter **-gff** should be a GFF2 file that in our test case includes 100 gene models used to "mock" train (80) geneid for the hymenoptera species _M.cingulata_ as well as to "mock" evaluate (20) the resulting parameter file (which should be named by the program as **"M.cingulata.geneid.optimized.param"**).  

The coordinates represented in this GFF2 file are genomic and correspond to the contigs and scaffolds in "M.cingulata.4training.fa". When training geneid for any species the user should provide GeneidTRAINer with a GFF2 with this format.  

**3. -fastas** M.cingulata.4training.fa  

The mandatory command line parameter **-fastas** should consist of a multi-FASTA file containing a few contigs/scaffolds of (in this test case) the hymenoptera species _M.cingulata_ which incorporate the 100 gene models used to train/evaluate geneid for this species.  When training geneid for any species the user should provide GeneidTRAINer with a multi-FASTA file corresponding to the GFF2 models selected under **-gff**  

**4. -results** ./output/  

The mandatory **-results** parameter tells the pipeline in which directory to store the results of the training process. In our test case that would be **./output/** (remember that this is relative to the working folder "/data" within the docker container and that in relation the the user's machine the output folder would be a sub-directory of `<userselecteddir>`)    


**5. -userdata** config.ext  _(optional)_

the optional command line option **-userdata** should have the path to a "user-configurable file" (_i.e._ config.ext) in which the user rather than the program selects a few of the parameters needed to generate an optimized geneid prediction file for your species of interest. 

Currently the user can select **minimum and maximum intron size**, **minimum and maximum intergenic distance** (values to be incorporated into the gene model portion of the parameter file). The user can also select the **start and end coordinates of the acceptor, donor, start (and for a few species branch) profiles** which are incorporated into the parameter file as PWMs (POSITION-WEIGHT MATRICES) of Markov models of order 1.  

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

**6. -reduced** no

The mandatory command line parameter which for our test case should be set to **no** for the first time the user runs the mock training analysis tells geneidTRAINer to run the entire pipeline. Once it has been run once for a particular species and training set the user can set the option to **yes** which will start the process at the point when geneidTRAINer generates PWMs or Markov 1 models of the splice sites and start codon. Setting **-reduced** to **yes** would only be useful if the user decides to change at least one of the gene model/splice site profile length parameters using the option **-userdata**  

**7. -v `<userselecteddir>`:/data -w /data**  (we recommend the user only change the `<userselecteddir>` and leave all else as is)  

the docker option **-v** mounts a user-selected directory/volume in a directory called **/data** within the the docker container and the option **-w** sets **/data** to be the working directory within the docker container.  

**8. -u $(id -u):$(id -g)**  

In a **queuing shared-file system** the option **-u** gives the docker container permissions to write in the user's directory system.  

########################################################################################

## IMPORTANT: An actual **training set** for geneidTRAINer should:  

a) be made up of at least **400-500** protein-coding gene models (and up to **~2500** sequences to keep the training process as short as possible)  

b) contain a large proportion of **multi-exonic** genes (in order for geneid to accurately model the splice sites)    

c) include only **non-overlapping** gene models (both on the same and opposite strands)   

d) contain only **complete** gene sequences (with a first, all internal and final exons and including a canonical start and stop codons in the case of multi-exon genes and canonical start and stop codons in the case of single-exon genes)    

e) be made up of sequences longer than at least **150-200** amino-acids  

f) be constituted by sequences previously aligned to a curated protein database (_i.e._ Uniprot90) using a program such as BLASTP to ensure that the sequences of the candidates correspond to actual protein-coding genes  

g) include sequences that overlap with the database proteins above over at least **90%** of their length     

#########################################################################################

The output files/directory of geneidTRAINer should be created in the path **selected by the user** (_i.e._ `<userselecteddir>`/output). These include several files that are in most cases not relevant to the user. The most important file is the geneid parameter file which can (in a full training protocol NOT this mock example) be used to predict sequences on your species of interest, in this case:

_**M.cingulata.geneid.optimized.param**_

The user can also find statistics on the training process by going the directory:

**`<userselecteddir>`/output/statistics_M.cingulata** (if the results dir is selected to be "./output") 

Importantly, the statistics file (_i.e._ **"1_Oct_Tue_10_13_training_statistics"**) includes a graphical ASCII representation of the nucleotide information content within the profiles for the start codon and splice sites which are obtained after the user runs the program the first time. By looking at the profiles the user can decide whether she wants to change their automatically selected start and end coordinates using the **"config.ext"** file on a subsequent execution (**"-reduced no"**).

The start and spice site profile logos representing the nucleotide information content around the start codon, donor and acceptor sites can be obtained from **`<userselecteddir>`/output/statistics_`<speciesname>`/plots_`<speciesname>`** (where $SPECIES=_M.cingulata_ in our test case):  

**Acceptor.pdf**  
**Donor.pdf**  
**Start.pdf**  

You will also be able to find a **gff2ps** (**M.cingulata.pdf** in our test case, otherwise `<speciesname>`.pdf) diagram representing all genes predicted in the evaluation scaffold built by geneidTRAINer.

## **additional information concerning running geneidTRAINer in the context of DOCKER**

The default command line for geneidTRAINer (no external config file) given the test files above is:

**docker run -u $(id -u):$(id -g) -v `<userselecteddir>`:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced no**

Where `<userselecteddir>` is the user-selected working directory which INSIDE the docker container is mounted as "/data". The command line above also assumes that the user created a directory called **"input"** in the `<userselecteddir>` where it placed the files used by "geneidTRAINer" (geneidtrainerdocker). The results are put into a directory called **"output"**  as selected by the user in this example (to be appended to the working directory `<userselecteddir>`).

the option **"-reduced no"** tells the program to run the training from the beginning. If after having trained geneid for the species of interest AT LEAST ONCE the user wishes to retrain it starting only at the point where the splice sites and start profile lengths are selected (based on each nucletide information content) it can do so by setting "-reduced" to YES (**-reduced yes**). 

_This will *ONLY* be useful when combined with using an external config file (i.e. -userdata .input/config.ext as described above) with user-selected profile start and end coordinates for any of the splice sites or startcodon (branch sites in a subser of fungi) and/or different minimum and maximum intron and intergenic sizes than those selected automatically by geneidTRAINer._    
However, the user can also set -reduced to NO (**-reduced no**) and still provide an external config file (**"-userdata ./input/config.ext**") with non-default values of minimim/maximum intron/intergenic sizes and start/end coordinates for splice sites/start profiles. 

Therefore if the user decides to use the **config file** the command line should instead be:

**docker run -u $(id -u):$(id -g) -v $PWD/:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced <no/yes> -userdata ./input/config.ext**

## IMPORTANT: IN ORDER TO ABORT THE TRAINING PROCESS WHILE IT IS RUNNING (WHICH WILL ALSO KILL THE DOCKER CONTAINER) THE USER CAN JUST PRESS CONTROL+C
