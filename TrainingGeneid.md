# TRAINING Geneid

# INTRODUCTION

Geneid is a widely used, well established, very fast and low carbon-footprint _ab initio_ gene prediction program used to find genes in anonymous genomic sequences designed to possess an hierarchical structure. In the first step, splice sites, (possibly) branch sites, start and stop codons are predicted and scored along the sequence using either Position Weight Arrays (PWAs) or Markov Models (of order 1 or 2) depending on the number of available sites. In a second step, exons are built from the sites. Exons are scored as the sum of the scores of the defining sites, plus the log-likelihood ratio of a Markov Model for coding DNA.  Finally, from the set of predicted exons, the gene structure is assembled, maximizing the sum of the scores of the assembled exons.

Geneid has been used as an _ab initio_ prediction tool in several International genome projects in the last several years. To give a few examples, geneid has been used in the annotation of mammalian species sequenced in the early to mid 2000s (Rat, Mouse, and Cow), of the model fish _Tetraodon nigroviridis_, of the chicken genome, of the model insect _Drosophila melanogaster_, and of the unicellular ciliate organism _Paramecium tetraurelia_. Geneid was also employed in the efforts to annotate a number of importan crop plant species over the years (_i.e_ tomato, wheat, melon and bean).  Moreover, our group has also collaborated with the Broad Institute (MIT) for years in the effort to annotate several fungal genomes of significance to Human health, and has an on-going collaboration with in the US and England involved in the annotation of several hymenoptera species (mostly social bees and wasps). It was also one of the programs used to anotate the genome of the endangered Iberian lynx and it is also being used to annotate an endangered bird species which can be found in Northern Russia and Southeast of Asia (_Eurynorhynchus pygmeus_). Geneid is also one of the _ab initio_ tools used in ENCODE research consortium, an International project to identify all functional elements in the human genome sequence and it has been used to help re-annotate an Alzheimer's disease-model rodent (_Octodon degus_). 

In order for geneid to predict genes optimally on a particular species it is generally necessary to build a parameter file specific for that species or related group of organisms.  To generate a species-specific matrix it is necessary to "train" the program and geneid parameter configurations already exist for a number of eukaryotic species (http://genome.crg.es/software/geneid/index.html#parameters).  Training basically consists of computing position weight matrices (PWMs) or Markov models for the splice (and sometimes branch) sites and start codons and deriving a model for coding DNA (generally a Markov model of order 4 or 5).  The resulting species-specific parameter file is subsequently optimized through the modification of two internal parameters that determine the final weight to be assigned to the exons used to build the genes (eWF) and ratio of exon weight to splice site statistics (oWF) respectively. Optimization is also achieved by determining the minimum and maximum intron length (and intergenic distances) that geneid will allow when predicting genes. 

The basic requirements for a training set are an annotation file (http://www.sanger.ac.uk/Software/formats/GFF/) and a set of FASTA-format sequences corresponding to the gene models in the annotation file.

Generally as few as 200 gene models could be used to build an accurate geneid parameter file, but generally a user would want to have a larger number of sequences (> 500 - >2500) to build an optimally accurate matrix and also to be able to set aside some of the gene models for testing purposes. The gene models should also include a large percentage of spliced genes so that statistics for the canonical splice sites can be accurately computed. The newly developed parameter file is evaluated against a set of sequences previously set aside for this purpose by the pipeline itself which randomly puts aside 20% of the input gene models for evaluation of the parameter file accuracy.

Until relatively recently most training of geneid for different species, and subsequent evaluation of the newly built parameter file required separately running a relatively large number of programs and scripts directly from a Unix command line, programs which were often “wrapped” in BASH shell scripting language. This training was done manually, mostly _in-house_, and involved over 30 AWK and PERL programs.   The training/evaluation process also required running directly, or through BASH shell scripting command calls, eight C-language programs (including geneid). This training strategy, while effective, was cumbersome and could require up several days to complete. It also demanded that the user be knowledgeable of the all the intricacies of the training process, and of the different programming languages being used, as script and parameter file modifications were almost always needed. 

In this document we describe the development of a PERL language integration tool (geneidTRAINer4docker.pl), which, in the context of a docker container, allows us to combine all the above-mentioned scripts/programs into a single pipeline-like script. While the original geneidTRAINer program versions were designed to be user-interactive at a few steps along the execution flow, for the purposes of this version of geneidTRAINer (within docker), user-intervention is not allowed except for for a few option that can be provided _a priori_ through a config file. Finally, the user is not required to have much knowledge of the training process itself.  The geneidTRAINer4docker.pl script must be run directly from a Unix-compatible command line in a machine containing the latest version of docker. 

# 1.0 REQUIREMENTS FOR BUILDING A GENEID TRAINING SET

The training set for geneid should:

a) be made up of at least **400-500** protein-coding gene models (and up to -optimally- **~2500** sequences). Adding more sequences is possible but will likely not results in improvements in the newly generated parameter file. The training will naturally take longer to complete the more sequences we start with.    

b) contain a large proportion of **multi-exonic** genes (in order for geneid to accurately model the splice sites)    

c) include only **non-overlapping** gene models (both on the same and opposite strands)   

d) contain only **complete** gene sequences (with a first, all internal exons and final exon, canonical start and stop codons in the case of multi-exon genes and canonical start and stop codons in the case of single-exon genes)    

e) be made up of sequences longer than at least **100-150** amino-acids  

f) be constituted by sequences previously aligned to a curated protein database (_i.e._ Uniprot90) using a program such as BLAST to ensure that the sequences of the candidates correspond to actual protein-coding genes  **(recommended)**  

g) include sequences that overlap with the database proteins above over at least **90%** of their length  **(recommended)**   

# 2.0 DESCRIPTION OF TRAINING SCRIPT (GeneidTRAINer4Docker.pl)

## 2.1 INPUT OPTIONS

The script has a number of input options:

_Usage: /scripts/geneidTRAINer4docker.pl -species `<speciesname>` -gff `<inputpath><gffname>` -fastas `<inputpath><fastasname>` -results `<results_dir>` -reduced `<yes/no>` -userdata `<inputpath><configfilename>`(optional) -branch `<inputpath><memeprofilefilename>[space]<memeprofilenumber>` (optional)_  

However, this script is executed in the context of a docker container:  

_docker run -u $(id -u):$(id -g) -v `<userselecteddir>`:/data -w /data geneidtrainerdocker -species `<speciesname>` -gff `<inputpath><gffname>` -fastas  `<inputpath><fastasname>`  -results `<results_dir>` -reduced  `<yes/no>` -userdata  `<inputpath><configfilename>` (optional) -branch `<inputpath><memeprofilefilename>[space]<memeprofilenumber>` (optional)_

The required options are:

#### 1) "-species": the name of the species being trained with the format: `<first letter of the genus name>`DOT`<species>`  (_i.e._ _H.sapiens_); 

#### 2) "-gff": the path to a GFF-format file (version2) containing the coordinates of the gene models to be used in the training process; an example of a GFF2 file which can be used to test the pipeline can be obtained from (https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2)

#### 3) "-fastas": a single (multi-)FASTA file with the DNA sequences of the gene models plus a good number of flanking nucleotides; an example of a FASTA file which can be used to test the pipeline can be obtained from (https://public-docs.crg.eu/rguigo/Data/fcamara/geneidtrainer/testing/M.cingulata.4training.fa).   

#### 4) "-results": The directory (`<results_dir>`) where the results should be stored, as a sub-directory of `<userselecteddir>` (user-defined working directory). (_i.e._ "./output").  

#### 5) "-reduced": (yes/no) an indication of whether the full pipeline should be run or, if the user has already trained for this specied and data set once,a “reduced” version executed from the point at which the user chooses the boundaries of the different splice site (and occasionally branch) profiles

The optional command line parameters are:

#### 6) "-userdata" (config file name and path): a "user-configurable" file in which the user rather than the program selects a few of the parameters needed to generate an optimized geneid prediction file for your species of interest. (see beginning of this document for more details) 

#### 7) "-branch": the path to, name of branch site profile and the number of the appropriate profile, if it has been previously produced (using the motif-finding program MEME) and should be used in the training process; (ONLY USED FOR TRAINING A SUB-SET OF SPECIES OF FUNGI EXPECTED TO POSSESS A CONSERVED BRANCH POINT) 

## 2.2 PARAMETER FILE BUILDING MODULES  

We have developed a few parameter file-building modules in PERL programming language that are intrinsic to the training pipeline. **Geneid::Param** allows the parameter file to be built on the fly as the different splice site and coding statistics are computed and optimizations completed. **Geneid::Isocore** is a handler for isocores in the geneid parameter file and is used mainly in the building of parameters for those species that possess more than one isocore. **Geneid::GeneModel** is a handler for GeneModel objects in GeneID parameter files. **Geneid::geneidCEGMA** contains a number of functions that we used to derive the coding potential of the candidate sequences using markov models of order 4 or 5. This module was based on same tool from HMMstar.pm (Ian Korf) and on a module used by the program CEGMA.

## 2.3 NEWLY DEVELOPED PARAMETER FILE AND TRAINING STATISTICS OUTPUT FILE  

The newly developed parameter file will be stored in the directory set by the mandatory **-results** input option. It will be named as follows: _**`<speciesname>`.geneid.optimized.param**_ (where the species name is taken from the mandatory input option **-species**; _i.e._ "_H.sapiens_")

As the newly developed parameter file for a particular species is being built a text file is written which includes detailed statistics on the training set as well as on the training process itself. It also indicates the name of the optimized parameter file and the accuracy performance estimation (produced by running geneid on the 20% of sequences set aside at the beginning of the training process). This file will be produced each time the pipeline is executed and stored under the directory set by the option **-results** (_i.e._ "./output" in the example in the README). The file will be stored in a directory called "statistics_`<speciesname>`" (_i.e._ statistics_M.cingulata). The statistics file name will include information with regard to its creation date and time (in hours and minutes): "_i.e._ 1_Oct_Tue_10_13_training_statistics". 

## 2.4 DIAGRAMS OF SPLICE SITES AND START CODON PROFILES AND GFF2PS PLOT OF THE ANNOTATION ON THE TEST SEQUENCES USING THE NEW PARAMETER FILE  

The start and spice site profile logos representing the nucleotide information content around the start codon, donor and acceptor sites can be obtained from **statistics_`<speciesname>`/plots_`<speciesname>`**:

**Acceptor.pdf**  
**Donor.pdf**  
**Start.pdf**  

You will also be able to find a **gff2ps** (http://genome.crg.es/software/gfftools/GFF2PS.html) plot (named **`<speciesname>`.pdf**) representing all genes predicted in the evaluation scaffold built by geneidTRAINer. 

## 3.0 CONCLUSIONS

The automated training pipeline program (in the context of a docker container) described here allows for a much faster and less cumbersome training of the gene predictor geneid than possible before. It also requires much less “know-how” on part of the user. Furthermore, in most instances the performance of the newly developed matrix is the same or very close to the accuracy obtained with the previous, slower, more user-intervention dependent training strategy.
