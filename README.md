# GeneidTRAINerDocker

Docker container that includes the perl pipeline that we use to automatically train our #in-house# _ab initio_ program geneid (http://genome.crg.es/software/geneid/). It also creates a set of files used to evaluate geneid and other gene prediction programs that used in our #in house# protein-coding genome annotation pipeline (which will automated using NextFlow in the course of the next few months)

---------------------------------------------------------------------------------------

First the user must download the project repository containing GeneidTRAINer as a docker container. If you have git installed in your system this can be done by typing in the command line:

git clone https://github.com/fcamara7/GeneidTRAINerDocker.git

## **Steps required for running geneidTRAINer in the context of a Docker container**

You must have a recent version of docker installed in your system (https://docs.docker.com/install/). 

Once docker is installed in you system and in order build the docker container GeneidTRAINerDocker type the following:

**docker build -t  geneidtrainerdocker .** (make sure **Dockerfile** is in the directory from which you are running this command)  

To obtain the actual command line for the training pipeline the user can type:

**docker run -it geneidtrainerdocker**  
_Usage: /scripts/geneidTRAINer4docker.pl -species `<speciesname>` -gff `<inputpath><gffname>` -fastas `<inputpath><fastasname>` -results `<results_dir>` -reduced `<yes/no>` -userdata `<inputpath><configfilename>`(optional) -branch `<inputpath><memeprofilefilename>[space]<memeprofilenumber>` (optional)_ 

### The "minimal" set of options needed to run geneidTRAINer in the context of docker are: 

**docker run -u $(id -u):$(id -g) -v `<userselecteddir>`:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced no -userdata ./output/config.ext (optional)**      

The example above uses **test data** found in ***https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing*** which contains a number of files that can be used to test the geneidTRAINer program contained within the distributed docker image, as well as a sample config file where the user can select some values that would override the automatic selections made by the GeneidTRAINer pipeline.

`<userselecteddir>` is the user-selected working directory which INSIDE the docker container is mounted as "/data". The sample command line above also assumes that the user creates a directory called **"input"** under `<userselecteddir>` where the files used by "geneidTRAINer" are placed (you may actually create a folder with any name). 

## Below we briefly describe the command line options used to run geneidTRAINer in the context of the _sample sequences_ provided as a test case:  

**1. -species** M.cingulata   

The mandatory command line parameter **-species** should have the name of the species being trained with the **first letter** of the Genus name and the **species** designation, with **one dot and no spaces** between genus letter and species name. In this test case that would be **"M.cingulata"**    

**2. -gff** M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2  

The mandatory command line parameter **-gff** should be a GFF2 file that in our test case includes 100 gene models used to "mock" train (80) geneid for the hymenoptera species _M.cingulata_ as well as to "mock" evaluate (20) the resulting parameter file (which should be named by the program as **"M.cingulata.geneid.optimized.param"**). In the example above this file was placed in the _user-created_ folder **./input** which is a sub-directory of `<userselecteddir>` as far as the user is concerned. The user can change the name of the input folder as she sees fit always taking in consideration that it will be under the `<userselecteddir>`.

The coordinates represented in this GFF2 file are genomic and correspond to the contigs and scaffolds in "M.cingulata.4training.fa". When training geneid for any species the user should provide GeneidTRAINer with a GFF2 with this format.  

**3. -fastas** M.cingulata.4training.fa  

The mandatory command line parameter **-fastas** should consist of a multi-FASTA file containing a few contigs/scaffolds of (in this test case) the hymenoptera species _M.cingulata_ which incorporate the 100 gene models used to train/evaluate geneid for this species.  When training geneid for any species the user should provide GeneidTRAINer with a multi-FASTA file corresponding to the GFF2 models selected under **-gff**. In the example above this file was placed in the _user-created_ folder **input** which is a sub-directory of `<userselecteddir>` as far as the user is concerned. The user can change the name of the input folder as she sees fit always taking in consideration that it will be under the `<userselecteddir>`.

**4. -results** ./output/  

The mandatory **-results** parameter tells the pipeline in which directory to store the results of the training process. In our test case that would be **./output/** (remember that this is relative to the working folder "/data" within the docker container and that as far as the user is concerned the output folder is a sub-directory of `<userselecteddir>`)    

**5. -reduced** no

The mandatory command line parameter which for our test case should be set to **no** for the first time the user runs the mock training analysis tells geneidTRAINer to run the entire pipeline. Once it has been **run at least once** for a particular species and training set the user can set the option to **yes** which will start the process at the point when geneidTRAINer generates PWMs or Markov 1 models of the splice sites and start codon. Setting **-reduced** to **yes** would only be useful if the user decides to change at least one of the gene model/splice site profile length parameters using the option **-userdata**  

**6. -userdata** config.ext  _(optional)_

the optional command line option **-userdata** should have the path to a "user-configurable file" (_i.e._ config.ext) in which the user rather than the program selects a few of the parameters needed to generate an optimized geneid parameter file the species being trained. 

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

**7. -v `<userselecteddir>`:/data -w /data**  (we recommend the user only change the `<userselecteddir>` and leave all else as is)  

the docker option **-v** mounts a user-selected directory/volume in a directory called **/data** within the the docker container and the option **-w** sets **/data** to be the working directory within the docker container.  

**8. -u $(id -u):$(id -g)**  

In a **shared-file cluster system** the option **-u** gives the docker container permissions to write in the user's directory system.  

## results produced by geneidTRAINer

The output files/directory of geneidTRAINer should be created in the path **selected by the user** (_i.e._ `<userselecteddir>`/`<results_dir>` - "output" in this case). These include several files that are generally not relevant to the regular user. The most important file is the optimized **geneid parameter** matrix which can (in a full training protocol NOT this mock example) be used to predict sequences on your species of interest. It is named by the pipeline as: 

_**`<speciesname>`.geneid.optimized.param**_  (using our test data, M.cingulata.geneid.optimized.param)  

The user can also find **statistics** files built during the training process by changing to the folder:

**`<userselecteddir>`/`<results_dir>`/statistics_`<speciesname>`**   

(where `<results_dir>` is "output" and  `<speciesname>`=_M.cingulata_ in this example)    

Importantly, the statistics file (_i.e._ **"1_Oct_Tue_10_13_training_statistics"**) includes a graphical ASCII representation of the nucleotide information content within the profiles for the start codon and splice sites which are obtained after the user runs the program the first time. By looking at the profiles the user can decide whether she wants to change their automatically selected start and end coordinates using the **"config.ext"** file on a subsequent execution (**"-reduced no"**).  

The start and spice site profile logos representing the nucleotide information content around the start codon, donor and acceptor sites can be obtained from the folder: 

**`<userselecteddir>`/`<results_dir>`/statistics_`<speciesname>`/plots_`<speciesname>`**  

(where `<speciesname>`=_M.cingulata_ and `<results_dir>` is "output" in our test case):

**Acceptor.pdf**  
**Donor.pdf**  
**Start.pdf**  

You will also be able to find a **gff2ps**-generated diagram (**M.cingulata.pdf** in our test case, otherwise `<speciesname>`.pdf) representing all genes predicted in the evaluation scaffold built by geneidTRAINer.

########################################################################################

## - IMPORTANT - An actual **training set** for geneidTRAINer should:  

a) be made up of at least **400-500** protein-coding gene models (and up to **~2500** sequences). Adding more sequences is possible but will most likely not results in improvements in the newly generated parameter file.     

b) contain a large proportion of **multi-exonic** genes (in order for geneid to accurately model the splice sites)    

c) include only **non-overlapping** gene models (both on the same and opposite strands)   

d) contain only **complete** gene sequences (with a first, all internal exons and final exon, canonical start and stop codons in the case of multi-exon genes and canonical start and stop codons in the case of single-exon genes)    

e) be made up of sequences longer than at least **100-150** amino-acids  

f) be constituted by sequences previously aligned to a curated protein database (_i.e._ Uniprot90) using a program such as BLAST to ensure that the sequences of the candidates correspond to actual protein-coding genes  **(recommended)**  

g) include sequences that overlap with the database proteins above over at least **90%** of their length  **(recommended)**  

#########################################################################################

## IMPORTANT: IN ORDER TO ABORT THE TRAINING PROCESS WHILE IT IS RUNNING (WHICH WILL ALSO KILL THE DOCKER CONTAINER) THE USER CAN JUST PRESS CONTROL+C
