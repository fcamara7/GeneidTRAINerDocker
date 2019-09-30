# GeneidTRAINerDocker

Docker container containing the perl program that we use to automatically train the ab initio program geneid. It also creates a set of files used to evaluate geneid and other gene prediction programs that used in our #in house# protein-coding genome annotation pipeline (which will automated using NextFlow in the next few months)  

To start the docker container GeneidTRAINerDocker type the following:

**docker build -t  geneidtrainerdocker .**

In order to see what the command line for the program type:

**docker run -it geneidtrainerdocker**

_Usage: /scripts/geneidTRAINer4docker.pl -species H.sapiens -gff <gffname> -fastas <fastasname> -results <results_dir> -reduced <yes/no> -userdata <configfilenamepath> (optional) -branch <pathmemefilename profile#> (optional)_

The minimal set of option to to run geneidTRAINer in the context of docker are: 

**docker run -u $(id -u):$(id -g) -v $PWD/:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced no**

The example above uses test data in https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing.

***https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing should contain a number of files that can be used to test the geneidtrainer program contained within the distributed docker image, as well as a sample config file where the user can select some values that would override the automatic selections by GeneidTRAINer:***

**1. M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2

a GFF2 file that includes 100 gene models used to "mock" train (80) geneid for the hymenoptera species M.cingulata as well as to "mock" evaluate (20) the resulting parameter file (which should be named by the program as "M.cingulata.geneid.optimized.param". The coordinates represented in this are genomic and correspond to the contigs and scaffolds in "M.cingulata.4training.fa".When training geneid for any species the user should provide GeneidTRAINer with a ggf2 with this format. 


**2. M.cingulata.4training.fa

File containing a few contigs/scaffolds of the hymenoptera species M.cingulata which incorporate the 100 gene models used to train/evaluate geneid for this species. 

**3. config.ext 

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

**M.cingulata.geneid.optimized.param

The user can also find statistics on the training process by going the directory

$PWD/output/statistics_M.cingulata (if the results dir is set to be "./output") 

Importantly, the statistics file (i.e. "22_51_14_27_8_119_5_269_0_training_statistics") should included the profiles for the start codon and splice sites which are obtained after the user runs the program the first time. By looking at the profiles the user can decide whether she/he wants to change their automatically selected start and end coordinates using the "config.ext" file.


####IMPORTANT###

***In order to run geneidTRAINer you must have docker installed on your machine***

The default command line for geneidTRAINer (no external config file) given the test files above is:

**docker run -u $(id -u):$(id -g) -v $PWD/:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced no 

Where $PWD is the user-selected working directory which inside the docker container is mounted as "/data". The command line above aslo assumes that the user created a directory called "input" in the $PWD where it placed the files used by "geneidTRAINer" (geneidtrainerdocker). The results are put into a directory called "output" also created by the user in the working directory $PWD

the option "-reduced no" tells the program to run the training from the beginning. If after having trained geneid for the species of interest AT LEAST ONCE the user wishes to retrain it starting only at the point where the splice site and start profile length is selected it can do so by setting "-reduced" to YES (-reduced yes). This will *ONLY* be useful when combined with using an external config file (i.e. -userdata .input/config.ext as described above) with user-selected profile start and end coordinates for any of the splice sites or startcodon (branch sites in a subser of fungi) and/or different minimum and maximum intron and intergenic sizes than those selected automatically by geneidTRAINer.    

the user can also set -reduced to NO (-reduced no) and provide an external config file (-userdata ./input/config with non-default values of minimim/maximum intron/intergenic sizes and start/end coordinates for splice sites/start profiles. 


If the user decides to use the config file the command line should instead be:

**docker run -u $(id -u):$(id -g) -v $PWD/:/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -reduced <no/yes> -userdata ./input/config.ext
#########################################################################


**The remainder of this document contains more detailed information concerning the gene prediction program geneid, the files required to run the geneid training pipeline (geneidTRAINer) in the context of the docker container, and the different input options.  


INTRODUCTION

Geneid is a widely used, well established, ab initio gene prediction program used to find genes in anonymous genomic sequences designed to possess an hierarchical structure. In the first step, splice sites, (possibly) branch sites, start and stop codons are predicted and scored along the sequence using either Position Weight Arrays (PWAs) or Markov Models (of order 1 or 2) depending on the number of available sites. In a second step, exons are built from the sites. Exons are scored as the sum of the scores of the defining sites, plus the log-likelihood ratio of a Markov Model for coding DNA.  Finally, from the set of predicted exons, the gene structure is assembled, maximizing the sum of the scores of the assembled exons.

Geneid has been used as an ab initio prediction tool in several International genome projects in the last several years. To give a few examples, geneid has been used in the annotation of mammalian species sequenced in the past few years (Rat, Mouse, and Cow), of the model fish T. nigroviridis, of the chicken genome, of the model insect Drosophila melanogaster, and of the unicellular ciliate organism Paramecium tetraurelia. Geneid was also employed in the efforts to annotate a number of importan crop plant species over the years (i.e tomato, melon and the bean).  Moreover, our group also collaborated with the Broad Institute (MIT) for years in the effort to annotate several fungal genomes of significance to Human health, and with a groups in the US and England involved in the annotation of several hymenoptera species (mostly social bees and wasps) . Geneid is also one of the ab initio tools used in ENCODE research consortium, an International project to identify all functional elements in the human genome sequence and it has been used to help re-annotate an Alzheimer's disease rodent (Octodon degus). 

In order for geneid to predict genes optimally on a particular species it is generally necessary to build a parameter file specific for that species or related group of organisms.  To generate a species-specific matrix it is necessary to "train" the program and geneid parameter configurations already exist for a number of eukaryotic species (link).  Training basically consists of computing position weight matrices (PWMs) or Markov models for the splice (and sometimes branch) sites and start codons and deriving a model for coding DNA (generally a Markov model of order 4 or 5).  The resulting species-specific parameter file is subsequently optimized through the modification of two internal parameters that determine the final weight to be assigned to the exons used to build the genes (eWF) and ratio of exon weight to splice site statistics (oWF) respectively. 

The basic requirements for a training set are an annotation file (http://www.sanger.ac.uk/Software/formats/GFF/) and a set of FASTA-format sequences corresponding to the gene models in the annotation file.

Generally as few as 100 gene models could be used to build a reasonably accurate geneid parameter file, but generally a user would want to have a larger number of sequences (> 500) to build an optimally accurate matrix and also to be able to set aside some of the gene models for testing purposes. The gene models should also include a large percentage of spliced genes so that statistics for the canonical splice sites can be accurately computed. The newly developed parameter file should be evaluated either against a set of sequences previously set aside for this purpose or, if not enough sequences were available, against the training set sequences themselves. When training and evaluation are perfomed on the same set of sequences a “ten-fold cross-validation” strategy is employed to reduce the performance bias that arises from training/evaluating on the same set of sequences. In this method 1/10 of the sequences are excluded from the training set, while the program is trained on the remaining sequences followed by the accuracy being estimated on the left-out sequences.

Until a few years ago most training of geneid for different species, and subsequent evaluation of the newly built parameter file required separately running a relatively large number of programs and scripts directly from a Unix command line, programs which were often “wrapped” in BASH shell scripting language. This training was done manually, mostly in-house, and involved over 30 AWK and PERL programs.   The training/evaluation process also required running directly, or through BASH shell scripting command calls, eight C-language programs (including geneid). This training strategy, while effective, was cumbersome and would require up to a week to complete. It also demanded that the user be knowledgeable of the all the intricacies of the training process, and of the different programming languages being used, as script and parameter file modifications were almost always needed. 

In this document we are going to describe the development of a PERL language integration tool (GeneidTrainer.pl), which allows us to combine all the above-mentioned scripts/programs into a single pipeline-like script. While the newly developed script was designed to be user-interactive at a few steps along the execution flow the user is not required to have much knowledge of the training process itself.  The GeneidTrainer.pl script must be run directly from a Unix command line. 

2. DESCRIPTION OF TRAINING SCRIPT (GeneidTRAINer.pl)

2.1 INPUT OPTIONS

The script has several input options. The required options are:

1) the name of the species being trained (i.e. “H.sapiens”); 

2) a GFF-format file (version2,link) containing the coordinates of the gene models to be used in the training process; an example of a GFF file which can be used to test the pipeline can be obtained from ().   

3) a single (multi-)FASTA file with the DNA sequences of the gene models plus a good number of flanking nucleotides; an example of a FASTA file which can be used to test the pipeline can be obtained from ().   

4) the name of the statistics file to be built as the the training pipeline is being executed; 

5) an indication of whether a branch site profile has been previously produced (using the motif-finding  program MEME) and should be used in the training process; (ONLY USED FOR SOME SPECIES OF FUNGI WHICH ARE EXPECTED TO POSSESS A CONSERVED BRANCH POINT) 

###6) an indication of whether the full pipeline should be run or, if the user has already trained once, or if a “reduced” version can be executed from the point at which the user chooses the boundaries of the different splice site (an occasionally branch) profiles; 

6) the user must also indicate the path to where she wants the results to be placed  


2.2 PARAMETER FILE BUILDING MODULES

We have developed two parameter file-building modules in PERL programming language that are intrinsic to the training pipeline. Geneid::Param allows the parameter file to be built on the fly as the different splice site and coding statistics are computed and optimizations completed. Geneid::Isocore is a handler for isocores in the geneid parameter file and is used mainly in the building of parameters for those species that possess more than one isocore.  

2.3 TRAINING STATISTICS OUTPUT FILE

As the newly developed parameter file for a particular species is being built a text file is written which includes detailed statistics on the training set as well as on the training process itself. It also indicates the name of the optimized parameter file and the accuracy  performance estimation.  


3.0 CONCLUSIONS

The semi-automated training pipeline script (in the context of a docker container) described in this paper allows for a much faster and less cumbersome training of the gene predictor geneid than it was possible before. It also requires much less “know-how” on part of the user. Furthermore, in most instances the performance of the newly developed matrix was the same or very close to the accuracy obtained with the previous, slower, more user-intervention dependent training strategy.
