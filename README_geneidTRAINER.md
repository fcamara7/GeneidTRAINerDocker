This document contains information concerning the gene prediction program geneid, the files required to run the geneid training pipeline (geneidTRAINer) in the context of a docker container, and all the different input options.  


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
