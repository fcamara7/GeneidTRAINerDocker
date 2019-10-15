#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Path;
use File::Basename;
use Geneid::Param;
use Geneid::Isocore;
use Geneid::geneid;
use Geneid::geneidCEGMA;
use Data::Dumper;

### This is the path to the awk scripts and C programs; change this to suit your system

# MAIN VARIABLES 
my $PROGRAM = "geneidTRAINer";
my $VERSION = "1.2docker";
my $path = "/scripts";
my $TMP = "./tmp/";
`mkdir -p $TMP`; 
my $TMPROOT   = "trainer_$$";
$SIG{TERM} = sub { die "Caught a sigterm $!" };

########PROGRAM SPECIFIC VARIABLES
my $count = 0;
my $gff = "";
my $fasta = "";
my $temptblcaps = "";
my $label = 0;
my $cutoff = -7; # for PWM profiles
my $pin = "-";
my $pout = "-";
my $results = "";
my $species = "";
my $answer = "";
my $validationanswer = "";
my $jacknifevalidate= "0";
my $useallseqs = 0;
my $tempgff4training = "";
my $tempgff4evaluation = "";
my $templocus_id = "";
my $templocus_id_new = "";
my $templocusid_eval = "";
my $seqsused = "";
my ($gpevalgff,$gpevalfa,$gpevaltbl) = ("","","");
my ($gpevalcontiggff,$gpevalcontigfa,$gpevalcontigtbl,$gpevalcontiglen) = ("","","","");
my ($outcdseval,$outintroneval,$outlocus_id_eval,$outgffeval,$inframeeval) = ("","","","",0);
my ($outcds,$outintron,$outlocus_id,$outgff,$inframe) = ("","","","",0);
my $total_seqs = "";
my $locus_id = "";
my $new_locus_id = "";
my $sortn = "sort -n";
my $seqs4training = "";
my $seqs4evaluation = "";
my $starttbl = "";
my $id;
my $tblseq = "";
my $usebranch = 0;
my @branchp = ();
my $memefile = "";
my $extconffile = "";
my $motifnumber = "";
my $memeanswer = "";
my ($branchmatrix,$prof_len_bra,$fxdbraoffset,$startbranch,$endbranch)=("","","","","");
my @evaluation = ();
my @jacknifeeval = ();
my ($bestIeWF,$bestIoWF,$bestAcc,$bestMin) = ("","","","");
my $fullengthbranchtbl = "";
my $fastasdir = "";
my $statsdir = "";
my $sitesdir = "";
my $introndir ="";
my $cdsdir = "";
my $geneidgffsorted = "";
my $total_genomic = "";
my $bckgrnd = "";
my $reduced = "";
my $reducedtraining = 0;
my ($outdonortbl,$totalnoncandon,$outacceptortbl,$totalnoncanacc,$outstarttbl,$totalnoncansta) =("","","","","","");
my ($gptraingff,$gptrainfa,$gptraintbl,$gptrainlen,$gpevallen) = ("","","","","");
my ($gptraincontiggff,$gptraincontigfa,$gptraincontigtbl,$gptraincontiglen) = ("","","","","");
my $totalseqs4training = "";
my $templist_train = "";
my $value = "";
my $gffseqseval = "";
my ($shortintron,$longintron,$minintergenic,$maxintergenic) = ("","","","");
my $donorsubprofile = "";
my $acceptorsubprofile = "";
my $startsubprofile = "";
my $branchsubprofile = "";
my $temp_jkf_geneid = "";
my $plotsdir = "";
my $contigopt = 0;
my $tempgeneidgffsorted = "";
my $tempgeneidgffsortedeval = "";
my $ext_flg = 0;
my $gfffiltered = "";
my $extdata = "";
my $useextdata = 0;
my $minintergenicusr = 0;
my $maxintergenicusr = 0;
my $shortintronusr = 0;
my $longintronusr = 0;
my $startusrdon = 0;
my $endusrdon = 0;
my $startusracc = 0;
my $endusracc = 0;
my $startusrsta = 0;
my $endusrsta = 0;
my $startusrbra = 0;
my $endusrbra = 0;
my $userprofile = 0;

GetOptions(
	   'species:s'         => \$species,
	   'gff:s'         => \$gff,
	   'fastas:s'      => \$fasta,
           'results:s'     => \$results,
           'reduced|red:s'	=> \$reduced,
           'userdata:s'      => \$extdata,
           'branch=s{2}'	=> \@branchp
	
	   	   );

my $usage = "Usage: $0 -species <speciesname> -gff <inputpath><gffname> -fastas <inputpath><fastasname> -results <results_dir> -reduced <yes/no> -userdata <inputpath><configfilename> (optional) -branch <inputpath><memeprofilefilename>[space]<memeprofilenumber> (optional)\n ";


print STDERR $usage and exit unless (($species =~ /^(\w\.\w+)|(\w+)$/i) && $gff && $fasta && $results && ($reduced =~ /^(yes|y)|(n|no)$/i)); 


####DECLARE A VARIABLE FOR A GIVEN SPECIES WHERE SEVERAL VALUES OBTAINED BY TRAINING NEED BE SAVED####


my $varsmemory = "${species}.variables";

print STDERR "\ndeclaring ${species}.variables file\n";

##########################################################################
##                                                                      ##
##                          INITIAL CHECKS                              ##
##                                                                      ##
##########################################################################

# Cheking if the external programs are in the path.
#C programs
system("which bash > /dev/null;")   && &go_to_die("bash is not found or is not executable");
system("which gawk > /dev/null;") && &go_to_die("gawk is not found or is not executable");
system("which egrep > /dev/null;") && &go_to_die("egrep is not found or is not executable");
system("which sed > /dev/null;") && &go_to_die("sed is not found or is not executable");
system("which geneid > /dev/null;")   && &go_to_die("geneid is not found or is not executable");
system("which SSgff > /dev/null;")  && &go_to_die("SSgff is not found or is not executable");
system("which shuf > /dev/null;") && &go_to_die("usort tool is not found or is not executable");
system("which pictogram > /dev/null;") && &go_to_die("The visualization program pictogram is not found or is not executable");
#BASH AWK
system("which gff2ps > /dev/null;")   && &go_to_die("The gff2ps package is not found or is not executable");
system("which gff2cds > /dev/null;")   && &go_to_die("The gff2cs package is not found or is not executable");
system("which gff2gp.awk > /dev/null;")   && &go_to_die("The gawk script gff2gp.awk is not found or is not executable");
system("which cds2gff > /dev/null;")   && &go_to_die("The gawk script cds2gff is not found or is not executable");
system("which frequency.awk > /dev/null;")   && &go_to_die("The gawk script frequency.awk is not found or is not executable");
system("which information.awk > /dev/null;")   && &go_to_die("The gawk script information.awk is not found or is not executable");
system("which submatrix.awk > /dev/null;")   && &go_to_die("The gawk script submatrix.awk is not found or is not executable");
system("which submatrix_order0.awk > /dev/null;")   && &go_to_die("The gawk script submatrix_order0.awk is not found or is not executable");
system("which Getkmatrix.awk > /dev/null;")   && &go_to_die("The gawk script Getkmatrix.awk is not found or is not executable");
system("which multiple_annot2one.awk > /dev/null;")   && &go_to_die("The gawk script multiple_annot2one.awk is not found or is not executable");
system("which logratio_kmatrix.awk > /dev/null;")   && &go_to_die("The gawk script logratio_kmatrix.awk is not found or is not executable");

system("which logratio_zero_order.awk > /dev/null;")   && &go_to_die("The gawk script logratio_zero_order.awk is not found or is not executable");
system("which preparedimatrixacceptor4parameter.awk > /dev/null;")   && &go_to_die("The gawk script preparedimatrixacceptor4parameter.awk is not found or is not executable");
system("which preparedimatrixdonor4parameter.awk > /dev/null;")   && &go_to_die("The gawk script preparedimatrixdonor4parameter.awk is not found or is not executable");
system("which preparetrimatrixstart4parameter.awk > /dev/null;")   && &go_to_die("The gawk script preparetrimatrixstart4parameter.awk is not found or is not executable");


########################################
########################################
###IF THERE IS A USER DATA CONFIG FILE##
########################################
########################################

if (-s $extdata) { ###IF THERE IS A USER DATA CONFIG FILE 

$useextdata = 1;
print "File exists and is named: ($extdata) \n\n";

open FILE, "<$extdata" || die "You need to provide a file with used derived data \n";
while(<FILE>){eval $_};
die "can't restore variables from $extdata: $@" if $@;
close FILE;


if ($startusrdon) {
 
    if (!$endusrdon || $endusrdon<=$startusrdon) {


    print STDERR "\nsince the start of the donor profile was set (\$startusrdon:$startusrdon) its end coordinate must be set (\$donorusrdon) with a value > than \$startusrdon:$startusrdon\nMPORTANT: The profile should be obtained automatically at least the first time so that the user can have access to the score distribution graph on the statistics output file\n" and exit;

}

}

if ($startusracc) {
 
    if (!$endusracc || $endusracc<=$startusracc) {

    print STDERR "\nsince the start of the acceptor profile was set (\$startusracc: $startusracc) its end coordinate must be set (\$endusracc) with a value > than \$startusracc:$startusracc\nIMPORTANT: The profile should be obtained automatically at least the first time so that the user can have access to the score distribution graph on the statistics output file\n" and exit;

}

}

if ($startusrsta) {
 
    if (!$endusrsta || $endusrsta<=$startusrsta) {

    print STDERR "\nsince the start of the start profile was set (\$startusrsta:$startusrsta) its end coordinate must be set (\$endusrsta) with a value > than \$startusrsta:$startusrsta\nIMPORTANT: The profile should be obtained automatically at least the first time so that the user can have access to the score distribution graph on the statistics output file\n" and exit;

}

}

if ($startusrbra) {
 
    if (!$endusrbra || $endusrbra<=$startusrbra) {

    print STDERR "\nsince the start of the branch profile was set (\$startusrbra:$startusrbra) its end coordinate must be set (\$endusrbra) with a value > than \$startusrbra:$startusrbra\nMPORTANT: The profile should be obtained automatically at least the first time so that the user can have access to the score distribution graph on the statistics output file\n" and exit;

}

}


if ($shortintronusr) {
 
    if (!$longintronusr || $longintronusr<=$shortintronusr) {

    print STDERR "\nsince a new minimum length for an intron was selected (\$shortintronusr:$shortintronusr) the longest allowed intron (\$longintronusr) must also be set with a value > than \$shortintronusr:$shortintronusr\n" and exit;

    }

}

if ($minintergenicusr) {
 
    if (!$maxintergenicusr || $maxintergenicusr<=$minintergenicusr) {

    print STDERR "\nsince a new minimum allowed intergenic distance was selected (\$minintergenicusr:$minintergenicusr) the longest intergenic distance (\$maxintergenicusr) must also be set with a value > than \$minintergenicusr:$minintergenicusr\n" and exit;

}

}

########################################
########################################
###IF THERE IS A MEME BRANCH PROFILE####
########################################
########################################
my $brel = 0;
$brel = @branchp;

if ($brel == "2") { ###IF THERE IS A MEME-DISCOVERED BRANCH POINT PROFILE

$usebranch = 1;

my $stringbranch = join(" ",@branchp);

$stringbranch =~ /^(.+?)\s+(\d)$/i ;

$memefile=$1;
$motifnumber=$2;

if (! -e "$memefile" || $motifnumber !~ /\d/i ) {
print "Please, indicate the name of the meme output file (make sure the file is in the right path) followed by the number of the motif the branch profile (i.e. write down -meme.txt 2- if the meme output file is called meme.txt and the second motif is the one corresponding to the branch profile) \n\n";

print STDERR $usage and exit;

}

}###IF THERE IS A MEME-DISCOVERED BRANCH POINT PROFILE END



if ($reduced =~ /^(no|n)/i)
{

    $reducedtraining = 0;
    print STDERR "should run the entire pipeline: $reducedtraining\n";
}

#############################################################
#############################################################
###REDUCED/SHORT TRAINING#####SKIPS ALL BUT BACKGROUND/PWM/MM5
#############################################################

if (-s "${results}$varsmemory" &&  $reduced =~ /^(yes|y)/i) { ###reduced/short training starting with PWMs 

$reducedtraining = 1;

print STDERR "\nshould run pipeline from PWM profile selection: $reducedtraining\n";

if (-d "$results") {

print STDERR "There is a directory named $results..\n\nbut elected the option to repeat the training for this species so NOT removing its contents\n";

}else{

    `mkdir -p $results;`;
;}


print STDERR "\n\nYou chose to continue the training of $species by assuming the cds, intronic sequences and splice sites have already been extracted \n\n";
  	     
open FILE, "<$results/$varsmemory" || die "You need to have run the training program once previously to execute the reduced version of the geneid training program \n";
while(<FILE>){eval $_};
die "can't restore variables from ${species}.variables: $@" if $@;
close FILE;


##CREATE A STATS FILE 
my @months = qw/Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec/;
my @days = qw/Sun Mon Tue Wed Thu Fri Sat Sun/;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my @timeData = ($mday,$months[$mon],$days[$wday],$hour,$min);

#STATS DIR WAS CREATED FIRST TIME PIPELINE IS RUN FOR A GIVEN SPECIES
my $statsout = $statsdir.join('_', @timeData)."_training_statistics";
###OPEN STATISTICS OUTPUT AT THIS TIME...EVERY TIME PIPELINE IS RUN
open SOUT,">$statsout";

if (!$useallseqs){print STDERR "\nThe reduced training process will use 80% of the gene-model sequences ($totalseqs4training)/20% will used for posterior evaluation of the newly developed parameter file ($gffseqseval)\n";} else {print STDERR "The reduced training process will use ALL of the gene-model sequences ($total_seqs)\n"};

print STDERR "\nA subset of $totalseqs4training sequences (randomly chosen from the $total_seqs gene models) was used for training\n";
print SOUT "GENE MODEL STATISTICS FOR $species\n\n";


} ########################## END OF SETTING SCRIPT FOR REDUCED TRAININING

##FULL TRAINING -MANDATORY FOR THE FIRST TIME GENEID IS TRAINED FOR A GIVEN SPECIES NOT REDUCED

if (!$reducedtraining) { #DO ONLY FIRST TIME YOU RUN FULL TRAINING PIPELINE 

if (-d "$results") {
print STDERR "There is a directory named $results..\nremoving directory and its contents\n";

	rmtree([ "$results/" ]);
print STDERR "\nmkdir -p $results\n";
        `mkdir -p $results;`;

}else{

    `mkdir -p $results;`;
;}

print STDERR "OPEN ${results}/$varsmemory\n";
open (STORV, ">${results}/$varsmemory")or die;

#######################################################

###MAKE A STATISTICS DIRECTORY

print STDERR "Create a statistics directory for this species\n";
`mkdir -p $results/statistics_${species}/`;
$statsdir = "$results/statistics_${species}/";

#####store statistics directory variable
print STORV Data::Dumper->Dump([$results], ['$results']);
print STORV Data::Dumper->Dump([$statsdir], ['$statsdir']);
####


##CREATE A STATS FILE 
my @months = qw/Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec/;
my @days = qw/Sun Mon Tue Wed Thu Fri Sat Sun/;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my @timeData = ($mday,$months[$mon],$days[$wday],$hour,$min);

#STATS DIR CREATED FIRST TIME PIPELINE IS RUN FOR A GIVEN SPECIES
my $statsout = $statsdir.join('_', @timeData)."_training_statistics";

###OPEN STATISTICS OUTPUT AT THIS TIME...AND EVERY TIME PIPELINE IS RUN
open SOUT,">$statsout";
#HEADER FOR STATS FILE
print SOUT "GENE MODEL STATISTICS FOR $species\n\n";

####Convert fasta to tabular format
print STDERR "\nConverting genomics fasta file ($fasta) to tabular format\n";

my $temptbl = ${results}.$species.".genomic.tbl";    

$temptbl = FastaToTbl($fasta,$temptbl);

     `sort -o $temptbl $temptbl`;
    
 print STDERR "\nactg to ACTG conversion of input fasta \n";

my $tblcaps = "";

  	    open LOCID, "gawk '{gsub(/_/,\"\",\$1);gsub(/\\./,\"\",\$1);print \$1, toupper(\$2)}' $temptbl |";
  	      while (<LOCID>) {
  		  $tblcaps .= $_;
  	      }
  	    close LOCID;

chomp $tblcaps;

$temptblcaps = ${results}.$species.".genomic.tbl";
            open FOUT, ">$temptblcaps";
 	    print FOUT "$tblcaps";
	    close FOUT;

        $value = `gawk '{print \$1}' $temptblcaps | sort | uniq | wc | gawk '{print \$1}'`; 
        chomp $value;

print STDERR "\nThe user has provided $value genomic sequences\n";

####store file with number of genomic sequences used in training
print STORV Data::Dumper->Dump([$value], ['$value']);
####store tabular file directory
print STORV Data::Dumper->Dump([$temptbl], ['$temptbl']);
####store CAPPED tabular  file directory
print STORV Data::Dumper->Dump([$temptblcaps], ['$temptblcaps']);
####

####CREATE FASTAS CDS; INTRON, SITES DIRs WITHIN PATH (ONLY FIRST TIME)
     print STDERR "\ncreate FASTAS, CDS, INTRON and SITES DIRECTORIES within PATH\n\n";
#CDS     
if (-d "$results/cds_${species}/") {
print STDERR "There is a directory named $results/cds_${species}/..\nRemove directory and its contents\n";
	rmtree([ "$results/cds_${species}/" ]);
        `mkdir -p $results/cds_${species}/;`;
        $cdsdir = "$results/cds_${species}/";

} else {
    `mkdir -p $results/cds_${species}/;`;
    $cdsdir = "$results/cds_${species}/";
}
#INTRON
if (-d "$results/intron_${species}/") {
print "There is a directory named $results/intron_${species}/!\nRemove directory and its contents\n";
	rmtree([ "$results/intron_${species}/" ]);
        `mkdir -p $results/intron_${species}/;`;
        $introndir = "$results/intron_${species}/";

} else {
    `mkdir -p $results/intron_${species}/;`;
    $introndir = "$results/intron_${species}/";
}
#SITES
   if (-d "$results/sites_${species}/") {
   print "There is a directory named $results/sites_${species}/!\nRemove directory and its contents\n";
   	rmtree([ "$results/sites_${species}/" ]);
           `mkdir -p $results/sites_${species}/;`;
           $sitesdir = "$results/sites_${species}/";

   } else {
       `mkdir -p $results/sites_${species}/;`;
       $sitesdir = "$results/sites_${species}/";
   }  
#FASTAS   
 if (-d "$results/fastas_${species}") {

print STDERR "There is already a directory named $results/fastas_${species}!\nRemove its contents and re-create it\n";

     rmtree([ "$results/fastas_$species/" ]);
    `mkdir -p $results/fastas_$species/`;
    $fastasdir = "$results/fastas_$species/";
    print STDERR "\n";
} else {
    print STDERR "Create a new one\nNo directory $results/fastas_${species} exists\n";
       `mkdir -p $results/fastas_$species/;`;
       $fastasdir = "$results/fastas_$species/";
}
#PLOTS JUST ONCE KEEP DATA IF IT HAD ALREADY BEEN CREATED
 if (-d "$statsdir/plots_${species}") { ###$statsdir = "$results/statistics_${species}/";

     print STDERR "\nThere is already a directory named $statsdir"."plots, however will keep contents, they may be overwritten";
    `mkdir -p $statsdir/plots_${species}`;
    $plotsdir = "$statsdir/plots_${species}";
    print STDERR "\n";
} else {
    print STDERR "Create a new one\nNo directory $statsdir"."plots_${species} exists\n";
       `mkdir -p $statsdir/plots_${species};`;
       $plotsdir = "$statsdir/plots_${species}";
}

#store fastas dir and plots dir
print STORV Data::Dumper->Dump([$fastasdir], ['$fastasdir']);
print STORV Data::Dumper->Dump([$plotsdir], ['$plotsdir']);
#######

       ##place genomic sequences in "fastas_$species" directory
       print STDERR "\nmove genomic sequences into \"fastas_$species\" directory\n";
       print STDERR "(also transfer genomic fasta length info)\n\n";
      ##do not create fastas in diretory if they are already created and their number corresponds to the number of sequences in thr array

###CONVERT GENOMICS FASTA TO MULTI FASTA AND PLACE THEM IN APPROPRIATE DIRECTORY

print STDERR "Convert $temptblcaps to multiple genomic fastas and place them in $fastasdir:\n";
       
       TblToFastaFile($fastasdir,$temptblcaps);
       print STDERR "\n\nConversion of $temptblcaps to multiple genomic fastas completed..\nAdd fasta sequence length information to same directory\n\n";

    opendir(D, $fastasdir) || die "Can't open directory: $!\n";
       
      my @files  = grep { 
            /^[\w+|\d+]/             # Begins with a character or number
	    && -f "$fastasdir/$_"   # and is a file
	} readdir(D);
  
       foreach my $file (@files) {
      
     
       	my $fnametable = ${file}.".tbl";
           $fnametable = FastaToTbl(${fastasdir}.$file,${fastasdir}.$fnametable);
	
          open FLEN, "<$fnametable";
              while (my $line = <FLEN>) {
 	     my ($name,$seq) = split (/\s+/,$line);
	       my $lengfasta = length($seq);
	      open FLEN2, ">$fastasdir/${name}_len";
                print FLEN2 "$name $lengfasta\n";
	      close FLEN2; 
           }
          close FLEN; 

   print STDERR "#";
	unlink 	$fnametable;	
       }
    closedir(D);

         print STDERR "\n";
#################################################

####get locus_id file only first time pipeline is run for a given species #ALL GENE MODELS 
    
     print STDERR "\nEliminate undesirable (_ and .) characters from $gff\n";
           
     my $filtergff = "";
     
     open LOCID, "gawk '{OFS=\"\\t\"}{gsub(/\\./,\"\",\$1);gsub(/\\./,\"\",\$9);gsub(/_/,\"\",\$0);print}' $gff |";
     while (<LOCID>) {
 	$filtergff .= $_;
     }
     close LOCID;

        
       open FOUT2, ">$gff";
       
       print FOUT2 "$filtergff";
       
       close FOUT2;
   
   print STDERR "\n $gff\n";
     
     
     print STDERR "\nObtain locus_id (list of genomic sequences / genes)\n";

     open LOCID2, "gawk '{print \$1,\$9}' $gff | sort | uniq |";
     while (<LOCID2>) {
 	$locus_id .= $_;
     }
     close LOCID2;
       
       $templocus_id = ${results}."/".$species."_locus_id";   
       
       open FOUT, ">$templocus_id";
       print FOUT "$locus_id";
       close FOUT;
    
###exit;     
  
########

###store locus id for a given species the first time we run the pipeline #ALL GENE MODELS
print STORV Data::Dumper->Dump([$templocus_id], ['$templocus_id']);
#####
 #####number of gene models TOTAL   
    $total_seqs = `gawk '{print \$2}' $templocus_id | sort | uniq | wc | gawk '{print \$2}'`;
    chomp $total_seqs;
 #####number of genomic sequences TOTAL
    $total_genomic = `gawk '{print \$1}' $templocus_id | sort | uniq | wc | gawk '{print \$1}'`;
    chomp $total_genomic;

    print STDERR "\nThe gff file ($gff) contains a total of $total_genomic genomic sequences and $total_seqs gene models\n";

##store total number of gene models and genomic sequences containing them
print STORV Data::Dumper->Dump([$total_seqs,$total_genomic], ['$total_seqs','$total_genomic']);

##### get a list of genes TOTAL
print STDERR "\nObtain list of all genes\n\n";
      my $list_seqs = "";
      open LOCID, "gawk '{print \$9}' $gff | sort | uniq |";
      while (<LOCID>) {
  	$list_seqs .= $_;
      }
      close LOCID;


  my $templist = ${results}."/".$species."_list_all_seqs";   
        open FOUT, ">$templist";
        print FOUT "$list_seqs";
        close FOUT;

#exit;  

####store list of gene models the first time the pipeline is run for a given species
print STORV Data::Dumper->Dump([$templist], ['$templist']);

} #FULL PIPELINE


###########################################################
###########################################################
#CREATING A PARAMETER FILE REGARDLESS OF WHETHER THE TRAINING IS COMPLETE OR SHORT VERSION
#########################################
########################################
#CREATE BLANK PARAMETER FILE############
my $param = Geneid::Param->new($species);
#set isochores to 1
$param->numIsocores(1);
$param->isocores([Geneid::Isocore->new()]);
########################################
########################################  
#######END CREATING A PARAMETER FILE REGARDLESS OF WHETHER THE TRAINING IS COMPLETE OR REDUCED
 
##############################################
 ####Select subset for training/evaluation###
##############################################

if (!$reducedtraining) { #RUN ONLY FIRST TIME FOR EACH SPECIES /ONLY FIRST TIME
      if ($total_seqs >= 25) { #switch to 400 in final distributed version	

 	   print "20% of the input sequences will be used for later evaluation of the accuracy of the parameter file\n";
	     
            $totalseqs4training = int(0.8*$total_seqs);
	  
  
	    print STORV Data::Dumper->Dump([$totalseqs4training], ['$totalseqs4training']);

	    print STDERR "\nA subset of $totalseqs4training sequences (randomly chosen from the $total_seqs gene models) was used for training\n";
	    #random pick (SHUF)
   	    open LOCID, "shuf $templocus_id | head -$totalseqs4training | sort | uniq |";
   	      while (<LOCID>) {
   		  $new_locus_id .= $_;
   	      }
   	    close LOCID;

  	    $templocus_id_new = ${results}."/".$species."_locus_id_training_setaside80";   
	    
	    open FOUT, ">$templocus_id_new";
  	    print FOUT "$new_locus_id";
  	    close FOUT;
	   
           print STORV Data::Dumper->Dump([$templocus_id_new], ['$templocus_id_new']);

 	   
###ASSUMING USER SELECTED TO SET ASIDE SEQUENCES FOR EVALUATION (20%)

	   $seqsused=`gawk '{print \$2}' $templocus_id_new | sort | uniq | wc | gawk '{ print \$1}'`;
           chomp $seqsused;

###################
####gff for training subset
####################
            my $gff4training = "";

            print STDERR "\nThe new training gff file includes $seqsused gene models (80% of total seqs)\n";
  	    open LOCID, "gawk '{print \$2\"\$\"}' $templocus_id_new | sort | uniq | egrep -wf - $gff |";
  	      while (<LOCID>) {
  		  $gff4training .= $_;          

  	      }
  	    close LOCID;

  	    $tempgff4training = ${results}."/".$species.".gff_training_setaside80";   
 	    open FOUT, ">$tempgff4training";
 	    print FOUT "$gff4training";
 	    close FOUT;
	    
	    print STORV Data::Dumper->Dump([$tempgff4training,$seqsused], ['$tempgff4training','$seqsused']);
	    

	    print STDERR "\nObtain list of training genes\n\n";
	    
	    my $list_seqs_train = "";
	    
	    open LOCID, "gawk '{print \$9}' $tempgff4training | sort | uniq |";
	    while (<LOCID>) {
		$list_seqs_train .= $_;
    }
     close LOCID;

	    $templist_train = ${results}."/".$species."_list_train_seqs_setaside80";   
	    open FOUT, ">$templist_train";
	    print FOUT "$list_seqs_train";
	    close FOUT;


#####Store variable with list of sequences set aside for training 
	   print STORV Data::Dumper->Dump([$templist_train], ['$templist_train']);
    

#########################  	    
####new locus_id for evaluation test set
#########################
   	   my $locusideval = "";
           open LOCID, "gawk '{print \$0\"\$\"}' $templist_train | egrep -vwf - $templocus_id |";
   	      while (<LOCID>) {      
   		  $locusideval .=$_;
   	      }
   	    close LOCID;
 	    chomp $locusideval;
 	    $templocusid_eval = ${results}."/".$species."_locus_id_evaluation_setaside20";
 	    open FOUT, ">$templocusid_eval";
 	    print FOUT "$locusideval";
 	    close FOUT; 

#####Store variable with list of sequences set aside for evaluating 
	    print STORV Data::Dumper->Dump([$templocusid_eval], ['$templocusid_eval']);

#########################
#####gff for evaluation test set
#########################

	    $gffseqseval=`gawk '{print \$2\"\$\"}' $templocusid_eval | sort | uniq | egrep -wf - $gff | gawk '{ print \$9}' | sort | uniq | wc | gawk '{print \$1}'`;
	    chomp $gffseqseval;

my $gff4evaluation = "";

            print STDERR "The evaluation gff file includes $gffseqseval gene models (20% of total seqs)\n\n";
  	    open LOCID, "gawk '{print \$2\"\$\"}' $templocusid_eval | sort | uniq | egrep -wf - $gff |";
  	      while (<LOCID>) {
  		  $gff4evaluation .= $_;
  	      }
  	    close LOCID;

  	    $tempgff4evaluation = ${results}."/".$species.".gff_evaluation_setaside20";   
 	    open FOUT, ">$tempgff4evaluation";
 	    print FOUT "$gff4evaluation";
 	    close FOUT;

###STORE INFO ON NUMBER OF SEQUENCES TO EVALUATE PLUS GFF FILE  OF SET ASIDE SEQUENCES..
	   print STORV Data::Dumper->Dump([$tempgff4evaluation,$gffseqseval], ['$tempgff4evaluation','$gffseqseval']);

$useallseqs = 0;
   	 
	    print STORV Data::Dumper->Dump([$useallseqs], ['$useallseqs']);
######################


    }# seqs > 400

  ####LOOP IF WE HAVE FEWER THAN 500 SEQUENCES
elsif ($total_seqs < 25) { # seqs < 400 FOR ALL INTENTS AND PURPOSES FOR THIS DOCKER VERSION OF THE PIPELINE JACKNIFFING IS DISABLED

    print STDERR "\nre-run pipeline making sure you have at least 400 sequences to train geneid with\n" and exit;



} # seqs < 25 (400)
    
}######### ABOVE ONLY EXECUTED FIRST TIME THE PIPELINE IS RUN FOR A GIVEN SPECIES



##############################################
################################CALL SUBS:####
##############################################


if (!$reducedtraining) { #ONLY FIRST TIME ("NOT SHORTER VERSION") FOR A GIVEN SPECIES

    if (!$useallseqs){ ##SET SEQS FOR EVAL AND TRAINING (SUBSETS)

####Convert general gff2 to geneid gff format
####extract and check cds and intron sequences. Remove inframe stops and check all seqs start with ATG and end with STOP

 ###TRAIN${results}.

     print STDERR "\nConvert general gff2 to geneid-gff format\n\n";  
       $tempgeneidgffsorted = generalGFFtoGFFgeneid($tempgff4training,$species,".train",$results);
        
    ($outcds,$outintron,$outlocus_id,$outgff,$inframe) = @{extractCDSINTRON($tempgeneidgffsorted,$templocus_id_new,".train",$results)};

   #  print STDERR " OUTSIDE EXTRACTCDSINTRON outgff: $outgff\noutlocus_id: $outlocus_id\n";
 ###TRAIN     


 ###EVAL
       $tempgeneidgffsortedeval = generalGFFtoGFFgeneid($tempgff4evaluation,$species,".eval",$results);
    
   ($outcdseval,$outintroneval,$outlocus_id_eval,$outgffeval,$inframeeval) = @{extractCDSINTRON($tempgeneidgffsortedeval,$templocusid_eval,".eval",$results)};
 ###EVAL

     #exit;
    	
} elsif ($useallseqs){ #USE SAME SEQS TO TRAIN/EVALUATE
####Convert general gff2 to geneid gff format
####extract and check cds and intron sequences. Remove inframe stops and check all seqs start with ATG and end with STOP
    
    print STDERR "\nConvert general gff2 to geneid-gff format $gff ###SAME SEQS USED TP TRAIN/EVALUATE\n\n";    
       
    $tempgeneidgffsorted = generalGFFtoGFFgeneid($gff,$species,".train");
    
    ($outcds,$outintron,$outlocus_id,$outgff,$inframe) = @{extractCDSINTRON($tempgeneidgffsorted,$templocus_id,".train",$results)};
      
} #USE SAME SEQS TO TRAIN/EVALUATE

####extract and check splice sites and start codon. Use only canonical info #IN SEQUENCES USED IN TRAINING
     ($outdonortbl,$totalnoncandon,$outacceptortbl,$totalnoncanacc,$outstarttbl,$totalnoncansta) = @{extractprocessSITES($outgff,$outlocus_id)};

#####prepare sequences for optimization of newly developed parameter file (TRAIN)    
    
    print STDERR "\nConvert gff to gp (golden-path-like)format (artificial contig - concatenated sequences - approx. 800 nt between sequences)\n";
     
    ($gptraincontiggff,$gptraincontigfa,$gptraincontigtbl,$gptraincontiglen) = @{processSequences4Optimization($outgff,".train",1,$results)};

    
    print STDERR "\nConvert gff to gp (golden-path-like)format (training set for later optimization -400-nt flanked sequences)\n";
     ($gptraingff,$gptrainfa,$gptraintbl,$gptrainlen) = @{processSequences4Optimization($outgff,".train",0,$results)};
    print STDERR "$gptraingff";

###STORE VARIABLE INFO IN DATA DUMPER###
    print STORV Data::Dumper->Dump([$outcds,$outintron,$outlocus_id,$outgff,$outcdseval,$outintroneval,$outlocus_id_eval,$outgffeval,$inframeeval,$tempgeneidgffsorted,$inframe,$outdonortbl,$totalnoncandon,$outacceptortbl,$totalnoncanacc,$outstarttbl,$totalnoncansta,$gptraingff,$gptrainfa,$gptraintbl,$gptrainlen,$gptraincontiggff,$gptraincontigfa,$gptraincontigtbl,$gptraincontiglen], ['$outcds','$outintron','$outlocus_id','$outgff','$outcdseval','$outintroneval','$outlocus_id_eval','$outgffeval','$inframeeval','$tempgeneidgffsorted','$inframe','$outdonortbl','$totalnoncandon','$outacceptortbl','$totalnoncanacc','$outstarttbl','$totalnoncansta','$gptraingff','$gptrainfa','$gptraintbl','$gptrainlen','$gptraincontiggff','$gptraincontigfa','$gptraincontigtbl','$gptraincontiglen']);
########################################


#NOT USING ALL SEQS FOR TRAINING/EVALUATION ALSO PROCESS EVAL SEQS
if (!$useallseqs){ 

#####prepare test set for evaluation of newly developed parameter file (EVAL)

    print STDERR "\nConvert gff to gp (golden-path-like)format (400-nt flanking)(test set for evaluation of new parameter file)\n";
     ($gpevalgff,$gpevalfa,$gpevaltbl,$gpevallen) = @{processSequences4Optimization($outgffeval,".eval",0,$results)};
    print STDERR "DONE\n";

    print STDERR "\nConvert gff to gp (golden-path-like)format (test set for evaluation of new parameter file - (artificial contig - concatenated sequences - approx. 800 nt between sequences)\n";
     ($gpevalcontiggff,$gpevalcontigfa,$gpevalcontigtbl,$gpevalcontiglen) = @{processSequences4Optimization($outgffeval,".eval",1,$results)};
print STDERR "DONE\n";

###STORE VARIABLE INFO IN DATA DUMPER###
   print STORV Data::Dumper->Dump([$gpevalgff,$gpevalfa,$gpevaltbl,$gpevallen,$gpevalcontiggff,$gpevalcontigfa,$gpevalcontigtbl,$gpevalcontiglen,$tempgeneidgffsortedeval], ['$gpevalgff','$gpevalfa','$gpevaltbl','$gpevallen','$gpevalcontiggff','$gpevalcontigfa','$gpevalcontigtbl','$gpevalcontiglen','$tempgeneidgffsortedeval']);
########################################

   
}#END OF NOT USING ALL SEQS




###DELETE DIRECTORIES NO LONGER NEEDED
print STDERR "the CDS-containing directory in $results/cds_${species}/ is no longer needed..\nRemove directory and its contents\n";
	rmtree([ "$results/cds_${species}/" ]);
print STDERR "the intron-containing directory in $results/intron_${species}/ is no longer needed..\nRemove directory and its contents\n";
	rmtree([ "$results/intron_${species}/" ]);

###EVERYTHING BELOW WILL BE RUN EVERYTIME THE TRAINING PIPELINE IS RUN WHETHER "REDUCED" OR "FULL" 

#####GET BACKGROUND SEQUENCES 
#User can set the length and number of the background sequences
my $kmer = 62;
my $numseqs = 100000;

             print "\nObtaining $numseqs background sequences of $kmer nucleotides each for estimating background frequencies of nucleotides\n";
	    
	    $bckgrnd = $species."_background.info";
	    $bckgrnd = getBackground($kmer,$fasta,$temptblcaps,$numseqs,$bckgrnd,$sitesdir);



  #  print STDERR "bckgrnd: $bckgrnd \n";

###STORE VARIABLE INFO IN DATA DUMPER###
   print STORV Data::Dumper->Dump([$sitesdir,$bckgrnd], ['$sitesdir','$bckgrnd']);
########################################
##############################################################
close STORV;

} ###IF NOT SHORT VERSION TRAINING TRAINING !$reducedtraining


####IF "REDUCED" TRAINING WILLS START AT THIS POINT: DONOR,ACCEPTOR,START AND BRANCH PROFILE STATS:



########
#########
#get donor site statistics
#########
#########
my $order = "0";
my $numbersites = `wc $outdonortbl | gawk '{print \$1}'`;
chomp $numbersites;
my $donoffset = "30"; #position before intron (last of exon (31) -1 for offset)

  if ($numbersites > 1400) {

   $order = "1";

  } elsif ($numbersites <= 1400){

   $order = "0";
}

print STDERR "\nThere are $numbersites donor sites, enough for a matrix of order $order, prior offset: $donoffset $outdonortbl $bckgrnd\n";

if ($startusrdon)

{$userprofile= 1;}

else {$userprofile = 0;}


my ($donormatrix,$prof_len_don,$fxddonoffset,$startdonor,$enddonor) = getKmatrix($outdonortbl,$bckgrnd,$order,$donoffset,1,0,0,0,0,0,0,$userprofile);

if (!defined @{$param->isocores}[0]->set_profile('Donor_profile',$prof_len_don,$fxddonoffset,$cutoff,$order,0,1,0,0,0,0,$donormatrix)){die "error in setting profile\n";}

my $donsub ="";

open LOCID, "gawk '{print  substr(\$2,($startdonor-3),($prof_len_don+6))}' $outdonortbl |";
  	      while (<LOCID>) {
  		  $donsub .= $_;
  	      }
  	    close LOCID;

  	    $donorsubprofile = $results.$species.".don.sub.profile";   
 	    open FOUT, ">$donorsubprofile";
 	    print FOUT "$donsub";
 	    close FOUT;

print STDERR "pictogram $donorsubprofile  $plotsdir/Donor -bits -land\n";

`pictogram $donorsubprofile $plotsdir/Donor -bits -land`;
`ps2pdf $plotsdir/Donor.ps $plotsdir/Donor.pdf`;
`rm $plotsdir/Donor.ps`;


########
#########
#get acceptor site statistics
#########
#########
 $order = "0";
 $numbersites = `wc $outacceptortbl | gawk '{print \$1}'`;
 chomp $numbersites;
 my $accoffset = "30"; #position after intron (first of exon (31) -1 for offset)

   if ($numbersites > 1400) {

    $order = "1";

   } elsif ($numbersites <= 1400){


    $order = "0";
 }

if ($startusracc)

{$userprofile= 1;}

else {$userprofile = 0;}


 print STDERR "\nThere are $numbersites acceptor sites, enough for a matrix of order $order, offset: $accoffset \n";


 my ($acceptormatrix,$prof_len_acc,$fxdaccoffset,$startacceptor,$endacceptor) = getKmatrix($outacceptortbl,$bckgrnd,$order,$accoffset,0,1,0,0,0,0,0,$userprofile);
if (!defined @{$param->isocores}[0]->set_profile('Acceptor_profile',$prof_len_acc,$fxdaccoffset,$cutoff,$order,0,1,0,0,0,0,$acceptormatrix)){die "error in setting profile\n";}


my $accsub ="";

open LOCID, "gawk '{print  substr(\$2,($startacceptor-3),($prof_len_acc+6))}' $outacceptortbl |";
  	      while (<LOCID>) {
  		  $accsub .= $_;
  	      }
  	    close LOCID;

  	    $acceptorsubprofile = $results.$species.".acc.sub.profile";   
 	    open FOUT, ">$acceptorsubprofile";
 	    print FOUT "$accsub";
 	    close FOUT;

print STDERR "pictogram $acceptorsubprofile  $plotsdir/Acceptor -bits -land\n";

`pictogram $acceptorsubprofile  $plotsdir/Acceptor -bits -land`;
`ps2pdf $plotsdir/Acceptor.ps $plotsdir/Acceptor.pdf`;
`rm $plotsdir/Acceptor.ps`;


#########
#########
#get start site statistics
#########
#########
 $order = "0";
 $numbersites = `wc $outstarttbl | gawk '{print \$1}'`;
 chomp $numbersites;
 my $staoffset = "30"; #before first position of the exon (31)minus 1 for offset)

   if ($numbersites > 5500) {

    $order = "2";

   } elsif ($numbersites <= 5500){

    $order = "0";
 }


if ($startusrsta)

{$userprofile= 1;}

else {$userprofile = 0;}


 print STDERR "\nThere are $numbersites start sites, enough for a matrix of order $order, offset: $staoffset \n";

 my ($startmatrix,$prof_len_sta,$fxdstaoffset,$startstart,$endstart) = getKmatrix($outstarttbl,$bckgrnd,$order,$staoffset,0,0,1,0,0,0,0,$userprofile);

####write to parameter file
if (!defined @{$param->isocores}[0]->set_profile('Start_profile',$prof_len_sta,$fxdstaoffset,$cutoff,$order,0,1,0,0,0,0,$startmatrix)){die "error in setting profile\n";}
#############################

my $stasub ="";

open LOCID, "gawk '{print  substr(\$2,($startstart-3),($prof_len_sta+6))}' $outstarttbl |";
  	      while (<LOCID>) {
  		  $stasub .= $_;
  	      }
  	    close LOCID;

  	    $startsubprofile = $results.$species.".sta.sub.profile";   
 	    open FOUT, ">$startsubprofile";
 	    print FOUT "$stasub";
 	    close FOUT;

print STDERR "pictogram $startsubprofile  $plotsdir/Start -bits -land\n";

`pictogram $startsubprofile  $plotsdir/Start -bits -land`;
`ps2pdf $plotsdir/Start.ps $plotsdir/Start.pdf`;
`rm $plotsdir/Start.ps`;


##OPTIONAL BRANCH STATS (FUNGI NORMALLY, AFTER RUNNING MEME)
if ($usebranch) {


$fullengthbranchtbl = processBranch($memefile,$motifnumber,$outintron);


########
#########
#get branch site statistics
#########
#########
 $order = "0";
 $numbersites = `wc -l $fullengthbranchtbl | gawk '{print \$1}'`;
 chomp $numbersites;
 my $braoffset = "32"; #before the A (branch) (33)minus 1 for offset)


if ($startusrbra)

{$userprofile= 1;}

else {$userprofile = 0;}

print STDERR "\nThere are $numbersites branch sites, enough for a matrix of order $order, offset: $braoffset \n";

my ($branchmatrix,$prof_len_bra,$fxdbraoffset,$startbranch,$endbranch) = getKmatrix($fullengthbranchtbl,$bckgrnd,$order,$braoffset,0,0,0,1,0,0,0,$userprofile);

####write to parameter file
if (!defined @{$param->isocores}[0]->set_profile('Branch_point_profile',$prof_len_bra,$fxdbraoffset,-50,$order,0,1,40,10,0,0,$branchmatrix)){die "error in setting profile\n";}
#############################

my $brasub ="";

open LOCID, "gawk '{print  substr(\$2,($startbranch-3),($prof_len_bra+6))}' $fullengthbranchtbl |";
  	      while (<LOCID>) {
  		  $brasub .= $_;
  	      }
  	    close LOCID;

  	    $branchsubprofile = $results.$species.".bra.sub.profile";   
 	    open FOUT, ">$branchsubprofile";
 	    print FOUT "$brasub";
 	    close FOUT;

print STDERR "pictogram $branchsubprofile  $plotsdir/Branch -bits -land\n";

`pictogram $branchsubprofile  $plotsdir/Branch -bits -land`;
`ps2pdf $plotsdir/Branch.ps $plotsdir/Branch.pdf`;
`rm $plotsdir/Branch.ps`;

} # USE usebranch


###DERIVE INITIAL/TRANSITION MARKOV MODEL

 my ($markovini,$markovtrans,$totalcoding,$totalnoncoding,$markovmodel) = @{deriveCodingPotential($outcds,$outintron)};

#add markov matrices to the parameter file 
if (!defined @{$param->isocores}[0]->Markov_order($markovmodel)){die "error in setting Markov_order\n";}
if (!defined @{$param->isocores}[0]->Markov_Initial_probability_matrix($markovini)){die "error in setting Markov_Initial_probability_matrix\n";}
if (!defined @{$param->isocores}[0]->Markov_Transition_probability_matrix($markovtrans)){die "error in setting Markov_Transition_probability_matrix\n";}
######################################
    

####PRODUCE FILE WITH STATS
if ($usebranch) {

($shortintron,$longintron,$minintergenic,$maxintergenic) = WriteStatsFile($species,$outintron,$outcds,$outgff,$inframe,$inframeeval,$seqsused,$totalnoncandon,$totalnoncanacc,$totalnoncansta,$markovmodel,$totalcoding,$totalnoncoding,$startdonor,$enddonor,$startacceptor,$endacceptor,$startstart,$endstart,$startbranch,$endbranch,$usebranch,$useallseqs);

} else {

($shortintron,$longintron,$minintergenic,$maxintergenic) = WriteStatsFile($species,$outintron,$outcds,$outgff,$inframe,$inframeeval,$seqsused,$totalnoncandon,$totalnoncanacc,$totalnoncansta,$markovmodel,$totalcoding,$totalnoncoding,$startdonor,$enddonor,$startacceptor,$endacceptor,$startstart,$endstart,0,0,0,$useallseqs);

} 



if (!$useextdata || ($shortintronusr==0 && $minintergenicusr==0) ) {

print STDERR "\nshortest intron: $shortintron\nlongest intron: $longintron\nminimum intergenic: $minintergenic\nmaximum intergenic: $maxintergenic\n";

} elsif ($useextdata && $shortintronusr>0 && $minintergenicusr>0) {

print STDERR "\nshortest intron: $shortintronusr\nlongest intron: $longintronusr\nminimum intergenic: $minintergenicusr\nmaximum intergenic: $maxintergenicusr\n";

} elsif ($useextdata && $shortintronusr>0 && $minintergenicusr==0) {

print STDERR "\nshortest intron: $shortintronusr\nlongest intron: $longintronusr\nminimum intergenic: $minintergenic\nmaximum intergenic: $maxintergenic\n";

} elsif ($useextdata && $shortintronusr==0 && $minintergenicusr>0) {

print STDERR "\nshortest intron: $shortintron\nlongest intron: $longintron\nminimum intergenic: $minintergenicusr\nmaximum intergenic: $maxintergenicusr\n";

} 


##################################################
###WRITE PRELIMINARY NON-OPTIMIZED PARAMETER FILE
$param->writeParam("$results/$species.geneid.param");
#print STDERR "PARAM: $param";
my $newparam = "$species.geneid.param";
print STDERR "newparam: $newparam";
################################################

################################
################################
###OPTIMIZE PARAMETER FILE######
################################
################################


print STDERR "\nOptimizing new parameter file\n\n";

my $opttype = "";

	   $contigopt= 1;

            print STDERR "\nThe parameter file will be optimized on an artificial contig made up of the concatenated flanked sequences (approx. 800 nt between genes)\n\n";
            print SOUT "\nThe parameter file will be optimized on an artificial contig made up of the concatenated flanked sequences (approx. 800 nt between genes)  \n\n";
	   $jacknifevalidate= 0;

####OPTIMIZATION FUNCTION NO BRANCH

my $array_ref = "";


##EXON WEIGHT PARAMETER 
     my $IeWF = "-6.5";
     my $deWF = "0.5";
     my $FeWF = "-2.0";
##EXON/OLIGO FACTOR PARAMETER     
     my $IoWF = "0.20";		
     my $doWF = "0.05";		
     my $FoWF = "0.70";
##Minimum Branch Profile Distance
     my $iMin = "5";		
     my $dMin = "2";		
     my $fMin = "11";
##ACCEPTOR CONTEXT
     my $iAccCtx = "40";		
     my $dAccCtx = "10";		
     my $fAccCtx = "70";

if (!$usebranch){ # no clear separate branch profile

###OPTIMIZATION FUNCTIONS

if (!$contigopt){

 @evaluation =@{OptimizeParameter($gptrainfa,$gptraingff,$newparam,0,0,0,0,$IeWF,$deWF,$FeWF,$IoWF,$doWF,$FoWF,0,0,0,0,0,0)};

($bestIeWF,$bestIoWF,$bestAcc,$bestMin,$array_ref) = @{BuildOptimizedParameterFile(\@evaluation,$usebranch,0,0,0)};

} elsif ($contigopt) {


@evaluation =@{OptimizeParameter($gptraincontigfa,$gptraincontiggff,$newparam,0,0,0,0,$IeWF,$deWF,$FeWF,$IoWF,$doWF,$FoWF,0,0,0,0,0,0)};

($bestIeWF,$bestIoWF,$bestAcc,$bestMin,$array_ref) = @{BuildOptimizedParameterFile(\@evaluation,$usebranch,0,0,0)};


}




} elsif ($usebranch){ #use separate branch profile

=head
my $respo = "";
do {
           print STDERR "Use automatically selected range values for the optimization of geneid eWF (exon weight)/oWF (exon/oligo factor)/Minimum Branch Distance from Acceptor and Context length in which Branch sites should be scored (AccCtx) internal parameters?\n\n(eWF: $IeWF to $FeWF; step $deWF\noWF: $IoWF to $FoWF; step $doWF\nMinBranchDistance: $iMin to $fMin; step $dMin\nAcceptorBranchCtx: $iAccCtx to $fAccCtx; step $dAccCtx)\n\nDo you prefer to change these values? (we do not recommed you change the Min Branch Distance and AccBranch context)";
  	   $respo = readline(STDIN);
       } while ($respo !~ /^(yes|y)|(n|no)$/i);

if ($respo =~/^(yes|y)/i) {
 	    
    my $sline = "";
    my $eline = "";
    my $dline = "";
            do {
            print STDERR "\nType new initial eWF (IeWF): ";
 	    $sline = readline(STDIN);
 	    } while ($sline !~/(-*[0-9]*\.*[0-9]+)/);
 	    $IeWF = $1;
 	    
	    do {
	    print STDERR "\nType new final eWF (FeWF): ";
 	    $eline = readline(STDIN);
 	    } while ($eline !~/(-*[0-9]*\.*[0-9]+)/|| $eline <= $sline);
 	    $FeWF = $1;
	    
	    do {
	    print STDERR "\nType step (delta) eWF (deWF)): ";
 	    $dline = readline(STDIN);
 	    } while ($dline !~/(-*[0-9]*\.*[0-9]+)/);
 	    $deWF = $1;
	    
	    do {
	    print STDERR "\nType new initial oWF (IoWF): ";
 	    $sline = readline(STDIN);
 	    } while ($sline !~/(-*[0-9]*\.*[0-9]+)/);
 	    $IoWF = $1;
 	    
	    do {
	    print STDERR "\nType new final oWF (FoWF): ";
 	    $eline = readline(STDIN);
 	    } while ($eline !~/(-*[0-9]*\.*[0-9]+)/ || $eline <= $sline);
 	    $FoWF = $1;
	    
	    do { 
	    print STDERR "\nType step (delta) oWF (doWF): ";
 	    $dline = readline(STDIN);
 	    } while ($dline !~/(-*[0-9]*\.*[0-9]+)/);
 	    $doWF = $1;
	    
	    do {
	    print STDERR "\nType new initial Min Branch Distance (iMin): ";
 	    $sline = readline(STDIN);
 	    } while ($sline !~/(-*[0-9]*\.*[0-9]+)/);
 	    $iMin = $1;
 	    
	    do {
	    print STDERR "\nType new final Min Branch Distance (fMin): ";
 	    $eline = readline(STDIN);
	    } while ($eline !~ /(-*[0-9]*\.*[0-9]+)/ || $eline <= $sline);
 	    $fMin = $1;
	    
	    do {
	    print STDERR "\nType step (delta) Min Branch Distance (dMin)): ";
 	    $dline = readline(STDIN);
 	    } while ($dline !~/(-*[0-9]*\.*[0-9]+)/);
 	    $dMin = $1;
	    
	    do {
	    print STDERR "\nType new initial Acceptor/Branch Context (iAccCtx): ";
 	    $sline = readline(STDIN);
	    } while ($sline !~/(-*[0-9]*\.*[0-9]+)/);
 	    $iAccCtx = $1;
 	    
	    do {
	    print STDERR "\nType new final Acceptor/Branch Context (fAccCtx): ";
 	    $eline = readline(STDIN);
	    } while ($eline !~ /(-*[0-9]*\.*[0-9]+)/ || $eline <= $sline);
	    $fAccCtx = $1;
	    
	    do {
	    print STDERR "\nType step (delta) Acceptor/Branch Context (dAccCtx): ";
 	    $dline = readline(STDIN);
 	    } while ($dline !~/(-*[0-9]*\.*[0-9]+)/);
 	    $dAccCtx = $1;

	}
=cut
####OPTIMIZATION FUNCTIONS

if (!$contigopt){

 @evaluation =@{OptimizeParameter($gptrainfa,$gptraingff,$newparam,1,$prof_len_bra,$fxdbraoffset,$branchmatrix,$IeWF,$deWF,$FeWF,$IoWF,$doWF,$FoWF,$iMin,$dMin,$fMin,$iAccCtx,$dAccCtx,$fAccCtx)};

($bestIeWF,$bestIoWF,$bestAcc,$bestMin,$array_ref) = @{BuildOptimizedParameterFile(\@evaluation,$usebranch,$prof_len_bra,$fxdbraoffset,$branchmatrix)};

}elsif ($contigopt) {

@evaluation =@{OptimizeParameter($gptraincontigfa,$gptraincontiggff,$newparam,1,$prof_len_bra,$fxdbraoffset,$branchmatrix,$IeWF,$deWF,$FeWF,$IoWF,$doWF,$FoWF,$iMin,$dMin,$fMin,$iAccCtx,$dAccCtx,$fAccCtx)};

($bestIeWF,$bestIoWF,$bestAcc,$bestMin,$array_ref) = @{BuildOptimizedParameterFile(\@evaluation,$usebranch,$prof_len_bra,$fxdbraoffset,$branchmatrix)};


}

#############################

}#USES BRANCH

my @evaluationinit = @$array_ref;
my @evaluationtest = ();

############
####EVALUATE PERFORMANCE OF NEW PARAMETER FILE ON TEST SET (IF PRESENT)
############ 

my $paramopt= "$species.geneid.optimized.param";

if (!$useallseqs){

#print STDERR "CHECK EVALUATE: $gpevalfa, $gpevalgff, $paramopt\n";

if (!$contigopt){

@evaluationtest =@{EvaluateParameter($gpevalfa,$gpevalgff,$paramopt)};

} elsif ($contigopt) {

@evaluationtest =@{EvaluateParameter($gpevalcontigfa,$gpevalcontiggff,$paramopt)};
}

if (!$usebranch) {

print STDERR "\nPerformance of new optimized parameter file on test set:\n\n".join("\t",@evaluationinit[2..$#evaluationinit]),"\n";
print SOUT "\nPerformance of new optimized parameter file on test set:\n\n".join("\t",@evaluationinit[2..$#evaluationinit]),"\n";
} elsif ($usebranch) {

print STDERR "\nPerformance of new optimized parameter file on test set:\n\n".join("\t",@evaluationinit[4..$#evaluationinit]),"\n";
print SOUT "\nPerformance of new optimized parameter file on test set:\n\n".join("\t",@evaluationinit[4..$#evaluationinit]),"\n";

}

print STDERR join("\t",@evaluationtest),"\n\n";
print SOUT join("\t",@evaluationtest),"\n\n";


} # if NOT using all seqs for training

if ($jacknifevalidate) {
###################################
###########SET UP FOR JACKNIFE: NOTE: NOT IMPLEMENTED IN THE DOCKER VERSION
###################################
#need output from processSequences4Optimization
 my $list4jkf = "";
# print STDERR "These are the sequences used for Jacknife (cross-validation) accuracy estimation\n";
  	    open LOCID, "grep '>' $gptrainfa| sed 's/>//g' | sort  |";
  	      while (<LOCID>) {
  		  $list4jkf .= $_;
  	      }
  	    close LOCID;

 my $templist4jkf = $species."_list_jacknife";   
 	    open FOUT, ">$templist4jkf";
 	    print FOUT "$list4jkf";
 	    close FOUT;

 my @list4jkf = split ("\n",$list4jkf);

    @list4jkf = sort @list4jkf;

my $totalseqstrain=`egrep '>' $gptrainfa -c | gawk '{print \$1}'`;
chomp $totalseqstrain;
######
#select how many sequences to group for jacknife
######
######SET UP FOR 10 FOLD VALIDATION
my $grpsiz = int (($totalseqstrain / 10) + 1);
######
#my $grpsiz = 1;
print STDERR "\nThe group size for 10 fold cross validation is $grpsiz\n\n";
print SOUT "\nThe group size for 10 fold cross validation is $grpsiz\n\n";
#deal with the fact that last array may have fewer elements than the grouping value
my $seqstofill = ($grpsiz - @list4jkf % $grpsiz);
if (@list4jkf % $grpsiz) {
  push @list4jkf, ("") x ($grpsiz - @list4jkf % $grpsiz)
}
#print STDERR "$outdonorbl"."\n";
#JACKNIFE FUNCTION

 if (!$usebranch) {

 @jacknifeeval = @{runJacknife (\@list4jkf,$seqstofill,$grpsiz,$outacceptortbl,$outdonortbl,$outstarttbl,$bckgrnd,$outcds,$outintron,$bestIoWF,$bestIeWF,$shortintron,$longintron,$minintergenic,$maxintergenic,$gptraintbl,$gptraingff,$startdonor,$enddonor,$startacceptor,$endacceptor,$startstart,$endstart,0,0,0,0,0,0)};

  } elsif ($usebranch) {

 @jacknifeeval = @{runJacknife (\@list4jkf,$seqstofill,$grpsiz,$outacceptortbl,$outdonortbl,$outstarttbl,$bckgrnd,$outcds,$outintron,$bestIoWF,$bestIeWF,$shortintron,$longintron,$minintergenic,$maxintergenic,$gptraintbl,$gptraingff,$startdonor,$enddonor,$startacceptor,$endacceptor,$startstart,$endstart,$usebranch,$fullengthbranchtbl,$startbranch,$endbranch,$bestAcc,$bestMin)};


 }


if (!$usebranch) {

print STDERR "\n\nPerformance of new parameter file after 10x cross-validation\n\n".join("\t",@evaluationinit),"\n";
print SOUT "\n\nPerformance of new parameter file after 10x cross-validation\n\n".join("\t",@evaluationinit),"\n";
} elsif ($usebranch) {

print STDERR "\n\nPerformance of new parameter file after 10x cross-validation\n\n".join("\t",@evaluationinit),"\n";
print SOUT "\n\nPerformance of new parameter file after 10x cross-validation\n\n".join("\t",@evaluationinit),"\n";

}

print STDERR join("\t",@jacknifeeval),"\n";
print SOUT join("\t",@jacknifeeval),"\n";

###############################################################
###############################################################

} # if jacknife

#########################################################################
######Obtain predictions on flanked sequences and plot them using gff2ps
#########################################################################


if ($jacknifevalidate && $useallseqs && !$contigopt) {

predictPlotgff2ps($paramopt,$gptrainfa,$gptraingff,$gptrainlen,$temp_jkf_geneid);

}elsif ($contigopt && !$useallseqs) {

predictPlotgff2ps($paramopt,$gpevalcontigfa,$gpevalcontiggff,$gpevalcontiglen,0);

}elsif ($contigopt && $useallseqs){


predictPlotgff2ps($paramopt,$gptraincontigfa,$gptraincontiggff,$gptraincontiglen,0);

}elsif (!$contigopt && !$useallseqs && $jacknifevalidate) {
   
predictPlotgff2ps($paramopt,$gptrainfa,$gptraingff,$gptrainlen,$temp_jkf_geneid);

}

print STDERR "\nplotting of the gff2ps graphs with the annotations and predictions given by the new $species parameter file completed\n";

#}
#####

#####unlink unnecessary files!


#####



close SOUT;

#remove tmpdir

rmtree([ "$TMP" ]);

######################################
#######################################
###END OF MAIN PORTION OF SCRIPT#######
#######################################
######BEGINNING OF SUBS################
#######################################
#######################################
#######################################
#######################################

#####FUNCTION TO EXTRACT CDS AND INTRONIC SEQUENCES AND GET RID OF TRANSLATIONS WITH IN-FRAME STOPS

sub extractCDSINTRON{

 	my ($gff,$locus_id,$type,$output) = @_;

# #####extract CDS and INTRON SEQUENCES

    	print STDERR "\nEXTRACT CDS and INTRON SEQUENCES from $type set..\n\n";
    	open(LOCUS,"<$locus_id");
     	my $count = 0;
    	while (<LOCUS>) {
    	    my ($genomic_id,$gene_id)=split;
	    `egrep -w '$gene_id\$' $gff > $TMP/$gene_id.gff`;
	    `SSgff -cE $results/fastas_$species/$genomic_id $TMP/$gene_id.gff | sed -e 's/:/_/' -e 's/ CDS//' >> $results/cds_${species}/${species}${type}.cds.fa`;
            `SSgff -iE $results/fastas_$species/$genomic_id $TMP/$gene_id.gff | sed -e 's/:/_/' -e 's/ Intron.*//' >> $results/intron_${species}/${species}${type}.intron.fa`;
   	    $count++;
   	    print STDERR "$count ..";
	   
    	}		       
    	close LOCUS;

    	print STDERR "DONE\n";

######tabulate CDS and INTRON SEQUENCES

   	print STDERR "\nCreate tabular format of CDS and INTRON sequences for $type sequences\n";
	

####CDS
	my $tempcdsfa = $cdsdir.${species}."$type".".cds.fa";
	print STDERR "$tempcdsfa\n\n";
        my $tempcds = ${output}.${species}."$type".".cds.tbl";
        $tempcds = FastaToTbl($tempcdsfa,$tempcds);
	print STDERR "cds tabular file created for $type sequences \n";
        
####INTRON
 	 my $tempintronfa = $introndir.${species}."$type".".intron.fa";
         my $tempintron = ${output}.${species}."$type".".intron.tbl"; 
         $tempintron = FastaToTbl($tempintronfa,$tempintron);

#####INTRONS LARGER THAN 0 ONLY

	    my $introntblpositive = "";
  	    open LOCID, "gawk '{if(length(\$2)>0){print \$1,\$2}}' $tempintron |";
  	    while (<LOCID>) {
	    
		$introntblpositive .= $_;
  	    }
	    close LOCID;
  	   
	    my $tempallintron_positive = ${output}.$species."$type".".intron_positivelength.tbl";   

	    open FOUT, ">$tempallintron_positive";
	    print FOUT "$introntblpositive";
	    close FOUT;
	
       print STDERR "intron tabular file created with introns with more than 0 nucleotides\n";

####GET LIST OF SEQUENCES WITH LENGTH >0 and EXCLUDE FROM CDS/locus_id/gff FILES SEQUENCES WITH INTRONS WITH 0 LENGTH
	
            my $intronzero = "";
  	    open LOCID, "gawk '{if(length(\$2)==0){print \$1}}' $tempintron | sed 's/\\(.*\\)\\..*/\\1\\_/' | sort | uniq |";
  	    while (<LOCID>) {
	    
		$intronzero .= $_;
  	    }
	    close LOCID;
  	   
	    my $tempall_intron_zero_list = ${output}.$species."$type".".intron_zerolength.list";   

	    open FOUT, ">$tempall_intron_zero_list";
	    print FOUT "$intronzero";
	    close FOUT;

	my $intronzero2 = "";

open LOCID, "gawk '{if(length(\$2)==0){print \$1}}' $tempintron | sed 's/\\(.*\\)\\..*/\\1/' | sort | uniq |";
  	    while (<LOCID>) {
	    
		$intronzero2 .= $_;
  	    }
	    close LOCID;
  	   
	    my $tempall_intron_zero_list2 = ${output}.$species."$type".".intron_zerolength.list2";   

	    open FOUT, ">$tempall_intron_zero_list2";
	    print FOUT "$intronzero2";
	    close FOUT;

#########FILTER SEQUENCES WITH 0 SIZE INTRONS FROM CDS!
	my $cdstblnozero = "";
  	    open LOCID, "egrep -vf $tempall_intron_zero_list $tempcds |";
  	    while (<LOCID>) {
		
		$cdstblnozero .= $_;
  	    }
	    close LOCID;
  	    
	    
	    my $tempallcds_nozero = ${output}.$species."$type".".cds_nozero.tbl";   

	    open FOUT, ">$tempallcds_nozero";
	    print FOUT "$cdstblnozero";
	    close FOUT;
###########ENSURE LOCUSID DOES NOT CONTAIN SEQUENCES WITH 0 SIZE INTRONS
	my $locusidnozero = "";
  	    open LOCID, "egrep -vwf $tempall_intron_zero_list2 $locus_id |";
  	    while (<LOCID>) {
		
		$locusidnozero .= $_;
  	    }
	    close LOCID;
  	    
	    
	    my $templocus_id_nozero= ${output}.$species."$type"."_locus_id_nozero";   

	    open FOUT, ">$templocus_id_nozero";
	    print FOUT "$locusidnozero";
	    close FOUT;
###########ENSURE GFF DOES NOT CONTAIN SEQUENCES WITH 0 SIZE INTRONS
	my $gffnozero = "";
  	    open LOCID, "egrep -vwf $tempall_intron_zero_list2 $gff |";
  	    while (<LOCID>) {
	    $gffnozero .= $_;
  	    }
	    close LOCID;
  	    
	    my $tempgffnozero = ${output}.$species."$type".".nozero.gff";   
	    
	    open FOUT, ">$tempgffnozero";
	    print FOUT "$gffnozero";
	    close FOUT;

  
####Convert sequences to protein format and check for in-frame stops
   	print STDERR "\nConvert sequences to protein format and check for in-frame stops and for proteins not starting with an M or not ending with a STOP\n\n";
       
        #SHOWS WHERE GENETIC CODE FILE IS LOCATED AND ITS NAME
        my $geneticcode = $path."/genetic.code";
	print STDERR "$path/genetic.code\n";
   	my $tempall_protein = ${output}.$species."$type".".protein";   
       	$tempall_protein = Translate($geneticcode,$tempallcds_nozero,$tempall_protein);
	
	my $inframestops=`gawk '{print \$2,\$1}' $tempall_protein | egrep '[A-Z]\\*[A-Z]\|^[^M]\|[^\\*] ' | gawk '{print \$2}' | wc | gawk '{print \$1}'`;
	chomp $inframestops;

  	print STDERR "\n\nWe found $inframestops sequences with in-frame stop signals/not starting with a methionine or not ending with a canonical stop codon \n\n";

  	###IF INFRAME######
  	if ($inframestops) {
  	    my $inframe = "";
	    my @inframe = ();
	    open LOCID, "gawk '{print \$2,\$1}' $tempall_protein | egrep '[A-Z]\\*[A-Z]|^[^M]|[^\\*] ' | gawk '{print \$2}' | sort | uniq |";
  	    while (<LOCID>) {
				
		push(@inframe, "$_" ); 	    

	    }
	    	   
	    close LOCID;
	   
	    foreach my $line (@inframe){
	        my (@frame)= split "_", $line;
		my $first = $frame[0];
	        $inframe .= "$first\n";
	    }
	    
	    my $inframe_protein = ${output}.$species."$type"."_INframe_NoMethionine_NoSTOP";
	    open FOUT, ">$inframe_protein";
	    print FOUT "$inframe";
	    close FOUT;
###REMOVE SEQUENCES WITH IN-FRAME STOPS FROM ORIGINAL CDS / INTRON / LOCUS_ID /GFF FILES AND PRINT NEW FILES
  	    print STDERR "\nremove sequences with in-frame stop signals from cds/intron_${species} files\n\n";
	    
	    
	    my $cdstbl2 = "";
  	    open LOCID, "sed 's/\\(.*\\)/\\1_/g' $inframe_protein | egrep -vf - $tempallcds_nozero |";
  	    while (<LOCID>) {
		
		$cdstbl2 .= $_;
  	    }
	    close LOCID;
  	    
	    
	    my $tempall_cds2 = ${output}.$species."$type".".cds_filter1.tbl";   

	    open FOUT, ">$tempall_cds2";
	    print FOUT "$cdstbl2";
	    close FOUT;
	    	    
	    my $introntbl2 = "";
  	    open LOCID, "sed 's/\\(.*\\)/\\1\.i/g' $inframe_protein | egrep -vf - $tempallintron_positive |";
  	    while (<LOCID>) {
	    
		$introntbl2 .= $_;
  	    }
	    close LOCID;
  	   
	    my $tempall_intron2 = ${output}.$species."$type".".intron_filter1.tbl";   

	    open FOUT, ">$tempall_intron2";
	    print FOUT "$introntbl2";
	    close FOUT;
	    
	    my $new_locus_id_filter1 = "";
	    
	    open LOCID, "sed 's/\\(.*\\)/\\1\$/g' $inframe_protein | egrep -vf - $templocus_id_nozero |";
  	    while (<LOCID>) {
	    $new_locus_id_filter1 .= $_;
  	    }
	    close LOCID;
  	    
	    my $templocus_id_new2 = ${output}.$species."$type"."_locus_id_filter_noinframe";   
	    
	    open FOUT, ">$templocus_id_new2";
	    print FOUT "$new_locus_id_filter1";
	    close FOUT;
	    
	    my $gffnew = "";

  	    open LOCID, "sed 's/\\(.*\\)_.*/\\1\$/g' $inframe_protein | egrep -vf - $tempgffnozero |";
  	    while (<LOCID>) {
	    $gffnew .= $_;
  	    }
	    close LOCID;
  	    
	    my $tempnewgff = ${output}.$species."$type".".noinframe.gff";   
	    
	    open FOUT, ">$tempnewgff";
	    print FOUT "$gffnew";
	    close FOUT;

######	    
	    return [$tempall_cds2,$tempall_intron2,$templocus_id_new2,$tempnewgff, $inframestops];


  	} else {		#######END IF THERE ARE INFRAME STOPS
  	    return [$tempallcds_nozero,$tempallintron_positive,$templocus_id_nozero,$tempgffnozero,0];
  	   }			#######END ELSE IF NO SEQS  ARE INFRAME


     }			#########sub extractCDSINTRON


####FUNCTION TO EXTRACT AND PROCESS SPLICE SITES AND START CODON

sub extractprocessSITES{

 	my ($gff,$locus_id) = @_;

 	#####SPLICE SITES
  	print STDERR "\nEXTRACT START AND SPLICE SITES\n\n";
 	my @newsites = ();
	open(LOC,"<$locus_id");
   	##print STDERR "$locus_id and $gff\n";
	
  	my $count = 0;
    	while (<LOC>) {
    	    my ($genomic_id,$gene_id)=split;
	  #  print STDERR "$genomic_id,$gene_id\n";
    	    `egrep -w '$gene_id\$' $gff > $TMP/$gene_id.gff`;
	  #  print STDERR "$gene_id $gff $TMP/$gene_id.gff \n\n";
    	    `SSgff -dabeE $results/fastas_$species/$genomic_id $TMP/$gene_id.gff > $TMP/${gene_id}.all_sites`;
    	    foreach my $site (qw(Acceptor Donor Stop Start)){
    	#	print STDERR "egrep -A 1 $site $TMP/${gene_id}.all_sites $sitesdir/${site}_sites.fa\n";
		`egrep -A 1 $site $TMP/${gene_id}.all_sites | sed -e '/--/d' -e '/^\$/d' >> $sitesdir/${site}_sites.fa`;
    	    }
    	$count++;
   	  print STDERR "site $count..";  
   	}			#while LOC
    	close LOC;
  	
	my $accsites="$sitesdir/Acceptor_sites.fa";
	#print STDERR "$accsites\n..";  
	my $donsites="$sitesdir/Donor_sites.fa";
	my $startsites="$sitesdir/Start_sites.fa";
	my $stopsites="$sitesdir/Stop_sites.fa";

	my $prestarttbl = "$sitesdir".$species.".canonical.start.tbl";
	$prestarttbl = FastaToTbl($startsites,$prestarttbl);
        my $acceptortbl = "$sitesdir".$species.".canonical.acceptor.tbl";
        $acceptortbl   = FastaToTbl($accsites,$acceptortbl);
	#print STDERR "$acceptortbl\n";
        my $donortbl =  "$sitesdir".$species.".canonical.donor.tbl";
        $donortbl = FastaToTbl($donsites,$donortbl);
      

##ADD N TO START SITES############

  	`gawk '{printf \$1" ";for (i=1;i<=60-length(\$2);i++) printf "n"; print \$2}' $prestarttbl > $sitesdir/${species}.canonical.start.complete.tbl`; 	 
  	my $starttbl = "$sitesdir".$species.".canonical.start.complete.tbl";
#################################
    
  	print STDERR "\n\nEliminate non-canonical donors/acceptors/starts:\n";

#  	##EXTRACT NON CANONICAL DONORS
  	my $noncanonical = "";
 	my $noncanonicalname = "";
  	my $totnoncanonical = "";
	my $totcanonical = "";
	my $newdonortbl = "";

  	open LOCID, "gawk '{print \$2}' $donortbl  | egrep -v '^[NATCGn]{31}GT' |";
  	 while (<LOCID>) { 
	  $noncanonical .= $_;
      }
  	close LOCID;
	
	my $tempdonornoncanonical = "$sitesdir".$species."_non_canonical_donor";   
 	    open FOUT, ">$tempdonornoncanonical";
 	    print FOUT "$noncanonical";
 	    close FOUT;

  	$totnoncanonical =`wc $tempdonornoncanonical | gawk '{print \$1}'`;
	chomp $totnoncanonical;
  	print STDERR "\nThere are $totnoncanonical non-canonical donors within the training set:\n";

  	###########################
  	if ($totnoncanonical) { #if there are non canonical donors
  	   
	    my @noncanonicalname = ();
	    open LOCID, "egrep -wf $tempdonornoncanonical $donortbl | gawk '{print \$1}' - | sort | uniq |";
  	     while (<LOCID>) {
		 
		 push(@noncanonicalname, "$_" );
		 
	}
  	    close LOCID;
	    
	    foreach my $line (@noncanonicalname){
	     #   my (@noncan)= split (/\.\d+:/, $line);
		my (@noncan)= split (/:/, $line);
		my $first = $noncan[0].":";
	        $noncanonicalname .= "$first\n";
		
	    }
          

	  #  unlink $tempdonornoncanonical;
    
	    my $tempnoncanonicalname = "$sitesdir".$species."_non_canonical_donor_seq_name";   
 	    open FOUT, ">$tempnoncanonicalname";
 	    print FOUT "$noncanonicalname";
 	    close FOUT;

  	    open LOCID, "egrep -vf $tempnoncanonicalname $donortbl |";
  	      while (<LOCID>) {
		  $newdonortbl .=$_;
	      }
  	    close LOCID;
	    
	    my $tempcanonicaldonor = "$sitesdir".$species.".canonical.donor.tbl";   
 	    open FOUT, ">$tempcanonicaldonor";
 	    print FOUT "$newdonortbl";
 	    close FOUT;
	    
	   # unlink $tempnoncanonicalname;
  	    

  	    my  $totcanonical =`wc $tempcanonicaldonor | gawk '{print \$1}'`;
	    chomp $totcanonical;

  	    print STDERR "\nThere are $totcanonical canonical donors within the training set:\n";
	    
	     push(@newsites, "$tempcanonicaldonor" );
	     push(@newsites, "$totnoncanonical" );
	}
	  
	else { #if there are no non-canonical
	    my  $totcanonical =`wc $donortbl | gawk '{print \$1}'`;
	    chomp $totcanonical;
  	    print STDERR "There are $totcanonical canonical donors within the training set:\n";
	    push(@newsites, "$donortbl" );
	    push(@newsites, "" );
  	    
  	}			#if there are no non-canonical
#  	###########################

#  	####
#  	##EXTRACT NON CANONICAL ACCEPTORS
  	 $noncanonical = "";
  	 $noncanonicalname = "";
	 $totcanonical = "";
	my $newacceptortbl = "";

  	open LOCID, "gawk '{print \$2}' $acceptortbl | egrep -v '^[NATCGn]{28}AG' |";
  	 while (<LOCID>) {
	     $noncanonical .= $_;
	 }
  	close LOCID;

	my $tempacceptornoncanonical = "$sitesdir".$species."_non_canonical_acceptor";   
 	    open FOUT, ">$tempacceptornoncanonical";
 	    print FOUT "$noncanonical";
 	    close FOUT;

  	$totnoncanonical =`wc $tempacceptornoncanonical | gawk '{print \$1}'`;
	chomp $totnoncanonical;
  	print STDERR "\nThere are $totnoncanonical non-canonical acceptors within the training set:\n";
  	###########################
  	if ($totnoncanonical) { #if there are non-canonical acceptors
  	    
	    my  @noncanonicalname = ();
	    open LOCID, "egrep -f $tempacceptornoncanonical $acceptortbl | gawk '{print \$1}' - | sort | uniq |";
	    while (<LOCID>) {
		 push(@noncanonicalname, "$_" );
	     }
		 close LOCID;

	    foreach my $line (@noncanonicalname){
	        my (@noncan)= split (/:/, $line);
		my $first = $noncan[0].":";
	        $noncanonicalname .= "$first\n";
		
	    }

	    unlink $tempacceptornoncanonical;


	    my $tempnoncanonicalname = "$sitesdir".$species."_non_canonical_acceptor_seq_name";   
 	    open FOUT, ">$tempnoncanonicalname";
 	    print FOUT "$noncanonicalname";
 	    close FOUT;

  	    open LOCID, "egrep -vf $tempnoncanonicalname $acceptortbl |";
  	      while (<LOCID>) {
		  $newacceptortbl .=$_;
	      }
  	    close LOCID;

	    #unlink $tempnoncanonicalname;
	    
	    my $tempcanonicalacceptor = "$sitesdir".$species.".canonical.acceptor.tbl";   
 	    open FOUT, ">$tempcanonicalacceptor";
 	    print FOUT "$newacceptortbl";
 	    close FOUT;
	    
	    #unlink $tempnoncanonicalname;
	    
  	    my  $totcanonical =`wc $tempcanonicalacceptor | gawk '{print \$1}'`;
	    chomp $totcanonical;

	    print STDERR "\nThere are $totcanonical canonical acceptors within the training set:\n";
	    
	    push(@newsites, "$tempcanonicalacceptor" );
	    push(@newsites, "$totnoncanonical" );

  	} 
	else {  #if there are only canonical use initial file list
	    my  $totcanonical =`wc $acceptortbl | gawk '{print \$1}'`;
	   # my  $totnoncanonical = "";
	    chomp $totcanonical;
  	    print STDERR "There are $totcanonical canonical acceptors within the training set:\n";
	    push(@newsites, "$acceptortbl" );
	    push(@newsites, "0" );
  	      	}			#if there are only canonical use initial file list
#  	###########################

#  	###
#  	##EXTRACT NON CANONICAL STARTS

  	$noncanonical = "";
  	$noncanonicalname = "";
	$totcanonical = "";
  	my $newstarttbl = "";

  	open LOCID, "gawk '{print \$2}' $starttbl | egrep -v '^[NATCGn]{30}ATG' |";
	while (<LOCID>) {
	    $noncanonical .= $_;
	}
  	close LOCID;

	my $tempstartnoncanonical = "$sitesdir".$species."_non_canonical_start";   
 	    open FOUT, ">$tempstartnoncanonical";
 	    print FOUT "$noncanonical";
 	    close FOUT;

	$totnoncanonical =`wc $tempstartnoncanonical | gawk '{print \$1}'`;
	chomp $totnoncanonical;

  	print STDERR "\nThere are $totnoncanonical non-canonical starts within the training set:\n";
  	###########################
  	
	if ($totnoncanonical) {	#if there are non-canonical starts
  	    
	    my @noncanonicalname = ();
	    open LOCID, "egrep -wf $tempstartnoncanonical $starttbl | gawk '{print \$1}' - | sort | uniq |";
  	      while (<LOCID>) {
		  push(@noncanonicalname, "$_" );
	  }
  	    close LOCID;

	    foreach my $line (@noncanonicalname){
	        my (@noncan)= split (/:/, $line);
		my $first = $noncan[0].":";
	        $noncanonicalname .= "$first\n";

	    }

	    unlink $tempstartnoncanonical;
	    
	    my $tempnoncanonicalname = "$sitesdir".$species."_non_canonical_start_seq_name";   
 	    open FOUT, ">$tempnoncanonicalname";
 	    print FOUT "$noncanonicalname";
 	    close FOUT;
	    
	    
  	    open LOCID, "egrep -vf $tempnoncanonicalname $starttbl |";
  	      while (<LOCID>) {
		  $newstarttbl .= $_;
	      }
  	    close LOCID;

	   # unlink $tempnoncanonicalname;
  	    
	    my $tempcanonicalstart = "$sitesdir".$species.".canonical.start.tbl";   
 	    open FOUT, ">$tempcanonicalstart";
 	    print FOUT "$newstarttbl";
 	    close FOUT;

	    #unlink $tempnoncanonicalname;

  	    my $totcanonical =`wc $tempcanonicalstart | gawk '{print \$1}'`;
	    chomp $totcanonical;

  	    print STDERR "\nThere are $totcanonical canonical starts within the training set:\n";
	    
	     push(@newsites, "$tempcanonicalstart" );
	     push(@newsites, "$totnoncanonical" );

	    
  	} 
	
	else {
	    my $totcanonical =`wc $starttbl | gawk '{print \$1}'`;
	    chomp $totcanonical;
  	    print STDERR "\nThere are $totcanonical canonical starts within the training set:\n";
	    push(@newsites, "$starttbl" );
	    push(@newsites, "0" );
	    
	}			#if there are only canonical starts
  	###########################

	return \@newsites;

     }				#subectractprocesssites


#######################################################3
###FUNCTION TO OBTAIN MARKOV MODELS CORRESPONDING TO THE CODING POTENTIAL
########################################################

 sub deriveCodingPotential{

     my $markov = "";
     my $markovm = "";
     my ($cds,$intron) = @_;

     my $totalcodingbases=`gawk '{ l=length(\$2); L+=l;} END{ print L;}' $cds`; 
     chomp $totalcodingbases;
     
     my $totalnoncodingbases=`gawk '{ l=length(\$2); L+=l;} END{ print L;}' $intron`;
     chomp $totalnoncodingbases;

     print STDERR "There are $totalcodingbases coding bases and $totalnoncodingbases non-coding bases on this training set:\n";
     
     if (($totalcodingbases>400000 && $totalnoncodingbases>100000) ||($totalcodingbases>375000 && $totalnoncodingbases>150000) ||($totalnoncodingbases>35000 && $totalcodingbases>(25*$totalnoncodingbases))) {
	 $markov = "5";
	 $markovm = "4";
	  print STDERR "Deriving a markov model of order $markov\n";

     } else {

	 $markov = "4";
	 $markovm = "3";
	 print STDERR "Deriving a markov model of order $markov\n";

     }

open (INIT,"<$intron");
       my @introndois =();
       while (<INIT>) {
	    my @i = split;
	    #print STDERR "SECOND FIELD $i[2]";
            push @introndois, $i[1];
	}
 	close INIT;

open (TRAN,"<$cds");
       my @coding =();
       while (<TRAN>) {
	    my @c = split;
            push @coding, $c[1];
	}
 	close TRAN;


     print STDERR "Intron model\n markov: ($markovm)";

my $intron_initial = new geneidCEGMA::SequenceModel('intron', 'FREQ', $markovm,
			 \@introndois, 10, 0);
 #print STDERR "($markovm) - intron initial: $intron_initial\n";
# $intron_initial->write("$DIR/intron.initial.5.freq");

my $intron_transition = new geneidCEGMA::SequenceModel('intron', 'MM', $markov,
			\@introndois, 10, 0);

# $intron_transition->write("$DIR/intron.transition.5.freq");

 
print STDERR "\nCoding model\n";

my $coding_initial = new geneidCEGMA::SequenceModel('coding', 'FREQ', $markov-1,
			 \@coding, 0.25, 2);

#$coding_initial->write("$path/coding.initial.5.freq");

my $coding_transition = new geneidCEGMA::SequenceModel('coding', 'MM', $markov,
			 \@coding, 0.25, 2);

#$coding_transition->write("$path/coding.transition.5.freq");


my $initial_logs = geneidCEGMA::log_ratio($coding_initial,$intron_initial);

my $transition_logs = geneidCEGMA::log_ratio($coding_transition,$intron_transition);

geneidCEGMA::write_log($initial_logs,"$results/coding.initial.5.logs.$species");

geneidCEGMA::write_log($transition_logs,"$results/coding.transition.5.logs.$species");
   #  print STDERR "\nINITIAL LOGS:".${$initial_logs}."\n";
#open (PROF1,"<$tempcdsintroninitgeneid");
open (PROF1,"<$results/coding.initial.5.logs.$species");
       my @profileinit =();
       while (<PROF1>) {
	    last if m/^\s/;
	    last if m/^[^ACGTacgt]/;
	    next if m/^#/;
	    chomp;
	    my @g = split;
            push @profileinit, \@g;
	}
 	close PROF1;

#open (PROF2,"<$tempcdsintrontrangeneid");
open (PROF2,"<$results/coding.transition.5.logs.$species");
       my @profiletran =();
       while (<PROF2>) {
	    last if m/^\s/;
	    last if m/^[^ACGTacgt]/;
	    next if m/^#/;
	    chomp;
	    my @g = split;
            push @profiletran, \@g;
	}
 	close PROF2;

return [\@profileinit,\@profiletran,$totalcodingbases,$totalnoncodingbases,$markov];

}#derive coding potential


###PROCESS SEQUENCES FUNCTION ( FLANKED GENE MODELS OBTAINED FOR OPTIMIZATION)

sub processSequences4Optimization{

my ($gff,$type,$contigopt,$output) = @_;

my $outtblname = "";
my $tblgp = "";
my $gff2gp = "";
my $fastagp = "";
my $gffgp = "";

open LOCID, "gff2gp.awk $gff | sort -k 1 |";
  	    while (<LOCID>) {
	    
		$gff2gp .= $_;
  	    }
	    close LOCID;

my $tempgff2gp = $output.$species.$type.".gp";   

	    open FOUT, ">$tempgff2gp";
	    print FOUT "$gff2gp";
	    close FOUT;
#print STDERR "BEFORE GETGENES: $fastasdir,$tempgff2gp,$path,$outtblname\n";
my $pretblgp = GetGenes($fastasdir,$tempgff2gp,$path,$outtblname);
#print STDERR "PRETBL AFTER GETGENES: $pretblgp \n";

print STDERR "\nGet sequences of 400-nt flanked sequences in tabular and gff formats\n";


  	    open LOCID, "gawk 'BEGIN{b=\"x\"}{if (\$1!=b){printf \"\\n\%s \",\$1}printf \"\%s\",\$6;b=\$1}END{printf \"\\n\"}' $pretblgp | sort | sed '/^\$/d' - | gawk '{print \$1, toupper(\$2)}' - |";
  	      while (<LOCID>) {
  		  $tblgp .= $_;
  	      }
  	    close LOCID;

my $tempgp_tbl = $output.$species.$type.".gp.tbl";   
 	    open FOUT, ">$tempgp_tbl";
 	    print FOUT "$tblgp";
 	    close FOUT;

            
             open LOCID, "gawk 'BEGIN{OFS=\"\\t\";pos=1;b=\"x\"}{if (\$1!=b){pos=1}; print \$1,\"annotations\",\$3,pos,pos+\$5-1,\"\.\",\"+\",\"\.\",\$1\$2; pos+=\$5;b=\$1 }' $pretblgp | egrep -v '(Intron|Utr)' - |";
  	      while (<LOCID>) {
  		  $gffgp .= $_;
  	      }
  	    close LOCID;


my $tempgp_gff = $output.$species.$type.".gp.gff";   
 	    open FOUT, ">$tempgp_gff";
 	    print FOUT "$gffgp";
 	    close FOUT;


print STDERR "DONE\n";

print STDERR "\nGet sequences of 400-nt flanked sequences in multi-fasta format\n";

my $tempgp_fa = $output.$species.$type.".gp.fa";   

$tempgp_fa = TblToFasta($tempgp_tbl,$tempgp_fa);

unlink $pretblgp;

print STDERR "\nSet up files for optimization\n\n";

my $seqslenggp = `gawk '{print \$1,length(\$2)}' $tempgp_tbl | sort -k1,1`;

my $tempseqlen = $output.$species.$type.".gp_cds_length";   
 	    open FOUT, ">$tempseqlen";
 	    print FOUT "$seqslenggp";
 	    close FOUT;
            

my $cdsgp = "";
           open LOCID, "gff2cds source=\"annotations\" $tempgp_gff | sort -k1,1 | join $tempseqlen - |";
  	      while (<LOCID>) {
  		  $cdsgp .= $_;
  	      }
  	    close LOCID;

my $tempcdsgp = $output.$species.$type.".cds_gp";   
 	    open FOUT, ">$tempcdsgp";
 	    print FOUT "$cdsgp";
 	    close FOUT;


my $gffgpeval = "";
           open LOCID, "gawk 'BEGIN{while (getline<ARGV[1]>0){len[\$1]=\$2;};ARGV[1]=\"\";OFS=\"\\t\";}{if (NR==1) {ant=\$1;print \$1,\$2,\"Sequence\",1,len[ant],\"\.\",\"\.\",\"\.\",\"\.\"};if (\$1!=ant) {print \"\#\$\";ant=\$1;print \$1,\$2,\"Sequence\",1,len[ant],\"\.\",\"\.\",\"\.\",\"\.\"}; print }' $tempcdsgp $tempgp_gff |";
  	      while (<LOCID>) {
  		  $gffgpeval .= $_;

  	      }
  	    close LOCID;

my $tempevalgpgff = $output.$species.$type.".gp_eval_gff";   
 	    open FOUT, ">$tempevalgpgff";
 	    print FOUT "$gffgpeval";
 	    close FOUT;


if ($contigopt){

#my $tempgp_fa = $species.$type.".gp.fa";   

#$tempgp_fa = FastaToTbl($tempgp_tbl,$tempgp_fa);


my @tabulargp = split(/\n/, $tblgp);
my $seq = "";
foreach my $line (@tabulargp){
	chomp $line;
	my @f = split " ", $line;
	$seq .= $f[1];
            }
my $lengp = length($seq);
my $foldedseqgp = fold4fasta($seq);
my $tempfastagpcontig = $output.$species.$type.".combined.gp.fa";
    open FOUT, ">$tempfastagpcontig";
    print FOUT ">$species\n$foldedseqgp\n";
    close FOUT;

my $temptabulargpcontig = $output.$species.$type.".combined.gp.tbl";
    open FOUT, ">$temptabulargpcontig";
    print FOUT "$species\t$seq\n";
    close FOUT;

my $seqslengcontiggp = `gawk '{print \$1,length(\$2)}' $temptabulargpcontig `;

my $tempseqlencontig = $output.$species.$type.".gp_cds_contig_length";   
 	    open FOUT, ">$tempseqlencontig";
 	    print FOUT "$seqslengcontiggp";
 	    close FOUT;


my $gpcontig = "";
open LOCID, "multiple_annot2one.awk species=$species leng=$lengp $tempcdsgp |";
  	    while (<LOCID>) {
	    
		$gpcontig .= $_;
  	    }
	    close LOCID;

my $tempgff2gpcontig = $output.$species.$type.".contig.gp.cds";   

	    open FOUT, ">$tempgff2gpcontig";
	    print FOUT "$gpcontig";
	    close FOUT;

my $cds2gffcontig = "";
 open LOCID, "cds2gff $tempgff2gpcontig | gawk 'BEGIN{OFS=\"\\t\";}{if (NR==1){print \"$species\",\"annotations\",\"Sequence\",\"1\",$lengp,\".\",\".\",\".\",\".\";print}else {print}}' -  | ";
         while (<LOCID>) {
  		  $cds2gffcontig .= $_;
  	      }
  	    close LOCID;

 my $tempgp_cdsgff_contig_eval= $output.$species.$type.".cds_gp_contig.eval.gff";   
 	    open FOUT, ">$tempgp_cdsgff_contig_eval";
 	    print FOUT "$cds2gffcontig";
 	    close FOUT;

return [$tempgp_cdsgff_contig_eval,$tempfastagpcontig,$temptabulargpcontig,$tempseqlencontig];

} elsif (!$contigopt)

{ 
    return [$tempevalgpgff,$tempgp_fa,$tempgp_tbl,$tempseqlen];
}

}#processSequences optimization



###GETGENES FUNCTION: EXTRACT FLANKED SEQUENCES FROM GENE MODELS FOR LATER OPTIMIZATION


sub GetGenes{

my ($path2gpath,$genesfname,$path,$outtblgp) = @_;
#print STDERR "IN FUNCTION: $path2gpath : $genesfname : $path : OUT: $outtblgp\n\n";

my $nonred=0;
my $onlynonred=0;
my $prevgene = "x";
my $prevchro = "x";
my $trail = "";
my %genenames;
$outtblgp="$results/outputtbl.txt";

chomp($path2gpath);
chomp($genesfname);
chomp($path);
chomp($outtblgp);

  if (!open(REFGENE,"< $genesfname")) {
  print "getgenes: impossible to open $genesfname\n";
  exit(1);
}

  if (!open(OUT,"> $outtblgp")) {
    print "getgenes: impossible to create $outtblgp\n";
    exit(1);
  }

while (<REFGENE>) {

  m/([\w\-\.:]+)\s+([\w\.\-:]+)\s+([\+\-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([^\n]+)/;

  my $name = $1;
  my $chro = $2;
  my $stra = $3;
  my $txSt = $4;
  my $txEn = $5;
  my $cdsS = $6;
  my $cdsE = $7;
  my $exoC = $8;
  my @exon = ($9 =~ m/(\d+)/g);
  my $cdsLe = $cdsE - $cdsS;
  my $txLe = $txEn - $txSt;
  my $cdsoffset = $cdsS - $txSt;
  my $redundant=0;
  my $i = 0; # exon counter
  my $j = 0;
  my $call = "";
  my $subseq = "";
  my $genomic = "";
  #my @tabular = ();

  if (!$onlynonred || ($onlynonred && !$redundant)) {
      open FLEN, "<$path2gpath${chro}_len";
      my $line = <FLEN>;
      chomp $line;
      my @le= split " ", $line;
      close FLEN;
###added code    
      my $chrotmptbl = "$TMP/tmp.tbl";
      $chrotmptbl = FastaToTbl($path2gpath.$chro,$chrotmptbl);
      #print STDERR "FATOTBL: $path2gpath"."$chro\n";
      open(IN,"<$chrotmptbl");
      my @tabular = ();
            while(<IN>){
	  #chomp; 
		#print STDERR "$_";
      push @tabular, "$_";
      }
      close IN;      
         #  print STDERR "\nGP: @tabular\n";
     # my @tabular = `FastaToTbl $path2gpath$chro`;
      my $subseq = "";
      #my $sublen = 0;
    foreach my $line (@tabular){
	chomp $line;
	my @f = split " ", $line;
	#print STDERR "$f[0]\n";
	$subseq .= $f[1];
        #$sublen += length($f[1]);
	#$countlines++;
	#if ($sublen >= ($numseqs * $k +1)){last;}
    }
######added

    
    if ($le[1] < $txEn) {
    
    my $newlen = $le[1];
    $genomic = substr($subseq, $txSt, ($newlen - $txSt));
    
    
}
    elsif ($le[1] >= $txEn) {
	$genomic = substr($subseq, $txSt, $txLe );

}

   # my $genomic = `$call`;
   # my $genomicLe = length($genomic);
      my $genomicLe = length($genomic);
      my $cdseq = "";

    if ($genomicLe == 0) {
      print STDERR "getgenes: gene of 0 length ($name), $call\n";
      next;
    }

  #  if ($genomicLe != $newlen) {
  #    print STDERR "getgenes: length mismatch ($name)\n";
  #    next;
  #  }



    for ($i=0;$i<$exoC;$i++) {

      my $utrB = 0;
      my $utrA = 0;
      my $utrS = 0;
      my $utrL = 0;
      my $exSt = $exon[$i] - $cdsS;
      my $exLe = $exon[$i+$exoC] - $exon[$i];
      my $exTy = "Internal";

      if ($exSt+$exLe > 0 && $exSt < $cdsLe) { # cds
  
        if ($exSt <= 0 || $i == 0) {
          if ($stra eq '+') {
            $exTy = "First";
          } else {
            $exTy = "Terminal";
          }
        }

        if ($exSt+$exLe >= $cdsLe || $i == $exoC-1) {
          if ($stra eq '+') {
            $exTy = "Terminal";
          } else {
            $exTy = "First";
          }
        }

        if ($exSt <= 0 && $exSt+$exLe >= $cdsLe) {
          $exTy = "Single";
        }

        if ($exSt < 0) {
          $utrB = 1;
          $utrS = $exSt;
          $utrL = abs($exSt);
          $exLe = $exLe - abs($exSt);
          $exSt = 0;
        }

        if ($exSt+$exLe > $cdsLe) {
          $utrA = 1;
          $utrS = $cdsLe;
          $utrL = $exLe - ($cdsLe - $exSt);
          $exLe = $cdsLe - $exSt;
        }

        my $iex;
        my $seq = substr($genomic,$exSt+$cdsoffset,$exLe);

        $seq = lc($seq);

        if ($stra eq '+') {  # forward

          if ($utrB) {
            my $iutr = $i+1;
            my $utrs = substr($genomic,$utrS+$cdsoffset,$utrL);
  
            $utrs = lc($utrs);
            $cdseq = $cdseq . "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n";
          }
  
          $iex = $i+1;
          $cdseq = $cdseq . "$name\t$chro\t$exTy\t$iex\t$exLe\t$seq\t$exon[$i]\t$exon[$i+$exoC]\n";

          if ($utrA) {
            my $iutr = $i+1;
            my $utrs = substr($genomic,$utrS+$cdsoffset,$utrL);
  
            $utrs = lc($utrs);
            $cdseq = $cdseq . "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n";
          }

        } else {             # reverse

          if ($utrB) {
            my $iutr = $exoC-$i;
            my $utrs = substr($genomic,$utrS+$cdsoffset,$utrL);
  
            $utrs = lc($utrs);
            $utrs =~ tr/acgt/tgca/;
            $utrs = reverse($utrs);
            $cdseq = "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n" . $cdseq;
          }

          $iex = $exoC-$i;
          $seq =~ tr/acgt/tgca/;
          $seq = reverse($seq);
          $cdseq = "$name\t$chro\t$exTy\t$iex\t$exLe\t$seq\t$exon[$i+$exoC]\t$exon[$i]\n" . $cdseq;

          if ($utrA) {
            my $iutr = $exoC-$i;
            my $utrs = substr($genomic,$utrS+$cdsoffset,$utrL);

            $utrs = lc($utrs);
            $utrs =~ tr/acgt/tgca/;
            $utrs = reverse($utrs);
            $cdseq = "$name\t$chro\tUtr\t$iutr\t$utrL\t$utrs\n" . $cdseq;
          }

        }

        if ($exTy ne "Single" && (($exTy ne "Terminal" && $stra eq '+') ||
           ($exTy ne "First" && $stra eq '-'))) {

          my $inSt = $exon[$i+$exoC] - $cdsS;
          my $inLe = $exon[$i+1] - $exon[$i+$exoC];

          if ($inSt+$inLe > 0 && $inSt < $cdsLe) {

            if ($inSt < 0) {
              print "getgenes: intron out of range! (1)\n";
              exit(1);
            }

            if ($inSt+$inLe > $cdsLe) {
              print "getgenes: intron out of range! (2)\n";
              exit(1);
            }

            $seq = substr($genomic,$inSt+$cdsoffset,$inLe);
            $seq = "\L$seq";

            my $iIn;

            if ($stra eq '+') { # forward
              $iIn = $j+1;
              $cdseq = $cdseq . "$name\t$chro\tIntron\t$iIn\t$inLe\t$seq\n";
            } else {
              $iIn = $exoC-$j-1;
              $seq =~ tr/acgt/tgca/;
              $seq = reverse($seq);
              $cdseq = "$name\t$chro\tIntron\t$iIn\t$inLe\t$seq\n" . $cdseq;
            }

          } else {
            print "getgenes.pl: intron out of range! (3)\n";
            if ($inSt+$inLe <= 0) {
              print "getgenes.pl: intron in 5' UTR\n";
            } else {
              print "getgenes.pl: intron in 3' UTR\n";
            }

            exit(1);
          }

          $j = $j + 1;
        }

      } else {  # UTRs
  
        $exSt = $exon[$i] - $txSt;
        $exLe = $exon[$i+$exoC] - $exon[$i];
  
        my $utrs = substr($genomic,$exSt,$exLe);
  
        if ($stra eq '+') {  # forward
          my $iutr = $i+1;

          $utrs = lc($utrs);
          $cdseq = $cdseq . "$name\t$chro\tUtr\t$iutr\t$exLe\t$utrs\n";
        } else {             # reverse
          my $iutr = $exoC-$i;

          $utrs = lc($utrs);
          $utrs =~ tr/acgt/tgca/;
          $utrs = reverse($utrs);
          $cdseq = "$name\t$chro\tUtr\t$iutr\t$exLe\t$utrs\n" . $cdseq;
        }

      }
    }

    print OUT $cdseq;

  } elsif ($onlynonred) {
    print STDERR "$name\n";
}

}
close(OUT);
return $outtblgp;

}#getgenes

##GET BACKGROUND SEQUENCES (ex. 62 Kmer) used to compute log likelihoods

sub getBackground {

my ($k,$fasta,$tbl,$numseqs,$background,$sited) = @_;

#my $fasta = $f;
#my $tbl = "";
my $countlines = 0;
my $totalseqs = `egrep -c \">\" $fasta`;
chomp $totalseqs;

open(IN,"<$tbl");
my @tabular = ();
while(<IN>){
     push @tabular, $_;
}
close IN;


print STDERR "\nThe total number of genomic sequences (in $fasta) is $totalseqs\n";

if ($totalseqs >= 1) {
    #print STDERR "in totalseqs if";
    my $seq = "";
    my $sublen = 0;
    foreach my $line (@tabular){
	chomp $line;
	my @f = split " ", $line;
	$seq .= $f[1];
        $sublen += length($f[1]);
	$countlines++;
	if ($sublen >= ($numseqs * $k +1)){last;}
    }
    #print STDERR "$seq";
    my $len = length($seq);
   # chomp $totalseqs;
    print STDERR "The length of the new fasta file is $len\n(concatenation of $countlines fastas (out of $totalseqs))\n";
    open FILE, ">${sited}$background";
    my $row = 1;
    print STDERR  "\nObtain background ($numseqs seqs) from fasta of $k nucleotides \n";
    for (my $n = 1;$n <= $numseqs;$n++) {
	my $kmer = "N";
	while ($kmer =~ m/[^ACGTacgt]/) {
	    my $rand = int (rand ($len - $k) + 1);
	    $kmer = substr($seq, $rand, $k);
	}
	#print STDERR  $n."..";
	print FILE $row++."\t$kmer\n";
    }
    print STDERR "FINISHED OBTAINING BACKGROUND INFO\n";
    close FILE;

}

return $sited.$background;

}# END getbackground function

sub processBranch {

my ($memefile,$motifnumber,$outintron) = @_;

#my $tempmotifbranch = $species."motif_branch";
#print STDERR "\nCHECK: begin to process branch (in function) meme: $memefile\n";
my $motif4branch = parseMEME($memefile,$motifnumber);
#print STDERR "\nCHECK: end of processing of branch (in function) meme: $memefile\n";
#open FOUT, ">$tempmotifbranch";
#	    print FOUT "$motif4branch";
#close FOUT;

my $realbranches = "";
     open LOCID, "gawk '{print \$4}' $motif4branch  |";
  	    while (<LOCID>) {
	       $realbranches .= $_;
  	    }
	    close LOCID;

my $temprealbranches = $species.".real_branches.tbl";

open FOUT, ">$temprealbranches";
	    print FOUT "$realbranches";
     close FOUT;


my $realbrancheslist = "";
     open LOCID, "gawk '{print \$1,\$4}' $motif4branch  | sed -e 's/^BL width.*//' -e /^\$/d |";
  	    while (<LOCID>) {
	       $realbrancheslist .= $_;
  	    }
	    close LOCID;


my $temprealbrancheslist = $species.".real_branches_list";

open FOUT, ">$temprealbrancheslist";
	    print FOUT "$realbrancheslist";
     close FOUT;

my $fullengthbranch = extractRealBranches($temprealbrancheslist,$outintron);

return $fullengthbranch;

}

sub parseMEME {

my ($memefi,$motifn) =  @_;
my $motiftemp = "Motif".$motifn;

if (-e "$motiftemp") {
			print STDERR "File $motiftemp exists already, erasing it..\n";
			unlink $motiftemp;
		}

open(MEMEFILE,"<$memefi")||die "Can't open $_: $!\n";
	my $in_blocks=0;
	while (<MEMEFILE>) {
		my $line=$_;
		open(MOTIF,">>$motiftemp")||die "Can't open $motiftemp: $!n";
		if (/^BL\s+MOTIF $motifn/) {
			$in_blocks=1;
		}
		if (/\/\// && $in_blocks) {
			$in_blocks=0;
			close MOTIF;
		}
		if ($in_blocks) {
			print MOTIF "$line";
		}
	    }

close MEMEFILE;

return $motiftemp;

}#END PARSE MEME

sub extractRealBranches {

my %branches;
my %sites;

my ($temprealbrancheslist,$outintron) = @_;

my $fullengthbranch = $species.".real.fullbranches.tbl";
my $fullengthbranchmod = "";

open(EXTRACT,">$fullengthbranch");
open(BRANCHES,"<$temprealbrancheslist")|| die"Can't open branches\n";

while (<BRANCHES>) {
	my ($id, $sequence)=split;
	$branches{$id}=$sequence;
}

close (BRANCHES);

open(INTRONS,"<$outintron")||die"Can't open introns\n";

while (<INTRONS>) {
	my ($id, $sequence)=split;
	$sequence=sprintf"%s%s%s",'N' x 30, $sequence, 'N' x 30;
	foreach my $key (%branches) {
		if (($id=~/$key/) && (my $pos=index ($sequence,$branches{$key}))) {
			$sites{$key}=substr($sequence,$pos-21,62);
			print EXTRACT "$id\t$sites{$key}\n";
		}
	}
}
close (BRANCHES);
close (EXTRACT);

open LOCID, "gawk '{if (\$2 !~ /^N/)print \$1,\$2}' $fullengthbranch |";
   	      while (<LOCID>) {      
   		  $fullengthbranchmod .=$_;
   	      }
   	    close LOCID;
 	    
 	    my $fullengthbranches = $species.".real.full_length_branches.tbl";
 	    open FOUT, ">$fullengthbranches";
 	    print FOUT "$fullengthbranchmod";
 	    close FOUT;


unlink $fullengthbranch;

return $fullengthbranches;


}



####GETKMATRIX FUNCTION (Splice site an Start codon PWMs)

sub getKmatrix {
 	my ($true_seqs,$false_seqs,$order,$offset,$donor,$accept,$star,$branch,$start,$end,$jacknife,$userprofile) = @_;
 	my $original_offset = $offset;
 	my @prof = ();
	my $tempinfolog;
  	my @orders = (qw(order-0 di tri order-4 order-5 order-6 order-7 order-8));
 	my $ordname = $orders[$order];
 	my $sort = "sort -n";
 	$sort = "sort -k1,1n -k2,2" if $order > 1;
 #	my @info = ($offset-1,$offset+1);
 	my $prof_len = 0;
 	my $info_thresh = ""; #bits
 	my $pid = $$;
 	my $true_seq_name = $true_seqs;
 	$true_seq_name =~ s/\.tbl$//;
 	#my $false_seq_name = $false_seqs;
 	#$false_seq_name =~ s/\.tbl$//;
 	my $false_seq_name = fileparse($false_seqs);
        #print STDERR "false_seq_name: $false_seq_name \n";
        ####Open true sequences
	#print STDERR "$true_seqs (true)\n";
 	open (T, "<$true_seqs");
 	$_ = <T>;
 	my @t = split;
 	my $len =  length($t[1]);
 	close T;
 	#######
 	####Open false sequences 
	#print STDERR "$false_seqs (false)\n";
 	open (F, "<$false_seqs") or die "Couldn't open $false_seqs: $!\n";
 	$_ = <F>;
 	my @f = split;
 	my $len2 =  length($f[1]);
 	close F;
 	######
 #	die "$len != $len2\n" if $len != $len2;
	#	print STDERR "$path/frequency.awk 1 $true_seqs > $true_seq_name.freq\n";
 	`gawk -f frequency.awk 1 $true_seqs > $true_seq_name.freq`;
	 #      print STDERR "$path/frequency.awk 1 $false_seqs > $sitesdir/$false_seq_name.freq\n";
 	`gawk -f frequency.awk 1 $false_seqs > $sitesdir/$false_seq_name.freq`;
	
	if ($donor) {
	`information.awk $sitesdir/$false_seq_name.freq $true_seq_name.freq | gawk 'NF==2 && \$1<=38 && \$1>=25' > $true_seq_name-$false_seq_name`;
	
	$tempinfolog = "$true_seq_name-$false_seq_name";
	print STDERR "tempinfolog: $tempinfolog \n";
	}
	if ($accept) {
	`information.awk  $sitesdir/$false_seq_name.freq $true_seq_name.freq | gawk 'NF==2 && \$1<=33 && \$1>=2' > $true_seq_name-$false_seq_name`;
         $tempinfolog = "$true_seq_name-$false_seq_name";
	}
	if ($star) {
	`information.awk $sitesdir/$false_seq_name.freq $true_seq_name.freq | gawk 'NF==2 && \$1<=37 && \$1>=25' > $true_seq_name-$false_seq_name`;
         $tempinfolog = "$true_seq_name-$false_seq_name";
	}
	
	if ($branch) {
	`gawk -f information.awk  $sitesdir/$false_seq_name.freq $true_seq_name.freq | gawk 'NF==2 && \$1<=41 && \$1>=28' > $true_seq_name-$false_seq_name`;
         $tempinfolog = "$true_seq_name-$false_seq_name";
	}
	

 	if (! $order) {
	     `gawk -f logratio_zero_order.awk $sitesdir/$false_seq_name.freq $true_seq_name.freq > $true_seq_name-log.$ordname-matrix`;		
	    

 	} else {
 	    `gawk -f Getkmatrix.awk $order $len $true_seqs | $sort > $true_seq_name.$ordname-matrix`;
 	    `gawk -f Getkmatrix.awk $order $len2 $false_seqs | $sort > $sitesdir/$false_seq_name.$ordname-matrix`;
 	    `gawk -f logratio_kmatrix.awk $sitesdir/$false_seq_name.$ordname-matrix $true_seq_name.$ordname-matrix > $true_seq_name-log.$ordname-matrix`;
 	}
 	#need to check output and then go on
 
	if (!$jacknife) {
	print STDERR "Information content profile\n";
    }
	
my $pwmstring = "";

	if ($donor && !$jacknife) {

$info_thresh = "0.15";
$pwmstring = "donor\n";

if (!$userprofile)
{


($start,$end) = BitScoreGraph($tempinfolog,$info_thresh,$offset,$pwmstring);



} else {

($start,$end) = BitScoreGraph($tempinfolog,$info_thresh,$offset,$pwmstring);

    $start = $startusrdon;
    $end = $endusrdon;
}

} elsif ($accept && !$jacknife) {

$info_thresh = "0.04";
$pwmstring = "acceptor\n";

if (!$userprofile)
{


($start,$end) = BitScoreGraph($tempinfolog,$info_thresh,$offset,$pwmstring);

} else {

($start,$end) = BitScoreGraph($tempinfolog,$info_thresh,$offset,$pwmstring);
    
    $start = $startusracc;
    $end = $endusracc;

}

} elsif ($star && !$jacknife) {

$info_thresh = "0.15";
$pwmstring = "start\n";

if (!$userprofile)
{
($start,$end) = BitScoreGraph($tempinfolog,$info_thresh,$offset,$pwmstring);
} else {

($start,$end) = BitScoreGraph($tempinfolog,$info_thresh,$offset,$pwmstring);

    $start = $startusrsta;
    $end = $endusrsta;

}


} elsif ($branch && !$jacknife) {

$info_thresh = "0.3";
$pwmstring = "branch\n";

if (!$userprofile)
{
($start,$end) = BitScoreGraph($tempinfolog,$info_thresh,$offset,$pwmstring);
} else {

($start,$end) = BitScoreGraph($tempinfolog,$info_thresh,$offset,$pwmstring);
    
    $start = $startusrbra;
    $end = $endusrbra;
}

}

###draw bit score bar graph function

sub BitScoreGraph {
	
        my ($infooutput,$info_thresh,$offset,$pwmstring) = @_;
	my @info = ($offset-1,$offset+1);
	print SOUT "$pwmstring";
        open (INFO,"<$infooutput");
 	while (<INFO>) {
 	    next if m/^#/;
 	    last if m/^\s/;
 	    last if m/^[^\d]/;
 	    chomp;
	    my @fields = split;
 	    printf STDERR "%2s %2.2f %s",($fields[0],$fields[1],"=" x int($fields[1] * 30)) ;
	    printf SOUT "%2s %2.2f %s\n",($fields[0],$fields[1],"=" x int($fields[1] * 30)) ;
 	    if ($fields[1] > $info_thresh) {
 		push (@info, $fields[0]) ;
 	    }
 	    print STDERR "\n";
 	}
 	close INFO;
 	my @sortedinfo = sort numerically @info;
 	my $start = (shift @sortedinfo);
 	$start = 1 if $start < 1;
 	my $end = pop @sortedinfo;
####


return ($start,$end);

    }#end BitScoreGraph

=head
if (!$jacknife){
	my $resp = "";
	my $sline = "";
	my $eline = "";
	do {
           print STDERR "#1 Automatically selected subsequence from $start to $end to use in profile.\nDo you prefer to change the start or end? ";
  	   $resp = readline(STDIN);
	  # exit;
       } while ($resp !~ /^(yes|y)|(n|no)$/i);
	if ($resp =~/^(yes|y)/i) {
 	    do {
	    print STDERR "\nType new start: ";
 	    $sline = readline(STDIN);
 	    } while ($sline !~/(\d+)/);
 	    $start = $1;
 	    do {
	    print STDERR "\nType new end: ";
 	    $eline = readline(STDIN);
 	    } while ($eline !~/(\d+)/ || $eline <= $sline);
 	    $end = $1;
 	}
    }#if not jacknife
=cut
 	$offset = $offset - $order;
 	#$start = $start - $order;
if (!$jacknife){
	$end = $end - $order;
}
 #	if ( $start < 1 ) {
 #	    $start = 1;
 #	}
 	$offset = $offset - $start + 1;
#print STDERR "end:$end offset:$offset start:$start\n";
#if (!$jacknife){
 #	print STDERR "new offset: $offset\nnew start: $start\nnew order: $order\n";
  #  }
	if ($order>=1 && $donor) {
	
	  my $preoffset = $offset + 2;
	  my $newoffset = $offset + 3;
	  my $posoffset = $offset + 4;
	    
 	`gawk -f submatrix.awk $start $end $true_seq_name-log.$ordname-matrix > $true_seq_name-log-info-pre.$ordname-matrix`;
	`preparedimatrixdonor4parameter.awk $preoffset $newoffset $posoffset $true_seq_name-log-info-pre.$ordname-matrix > $true_seq_name-log-info.$ordname-matrix `;
     # print STDERR "submatrix.awk $start $end $true_seq_name-log.$ordname-matrix/$true_seq_name-log-info.$ordname-matrix";

      } elsif ($order>=1 && $accept) {

          my $preoffset = $offset - 1;
	  my $newoffset = $offset;
	  my $posoffset = $offset + 1;
	    
 	`gawk -f submatrix.awk $start $end $true_seq_name-log.$ordname-matrix > $true_seq_name-log-info-pre.$ordname-matrix`;
	`preparedimatrixacceptor4parameter.awk $preoffset $newoffset $posoffset $true_seq_name-log-info-pre.$ordname-matrix > $true_seq_name-log-info.$ordname-matrix `;
    #	  print STDERR "submatrix.awk $start $end $true_seq_name-log.$ordname-matrix/$true_seq_name-log-info.$ordname-matrix";
 
      }  elsif ($order>=2 && $star) {

          my $preoffset = $offset - 2;
	  my $newoffset = $offset - 1;
	  my $posoffset = $offset;
	    
 	`gawk -f submatrix.awk $start $end $true_seq_name-log.$ordname-matrix > $true_seq_name-log-info-pre.$ordname-matrix`;
	`preparetrimatrixstart4parameter.awk $preoffset $newoffset $posoffset $true_seq_name-log-info-pre.$ordname-matrix > $true_seq_name-log-info.$ordname-matrix `; 
}
else {

   # print STDERR "submatrix_order0.awk $start $end $true_seq_name-log.$ordname-matrix\n";
	
	`gawk -f submatrix_order0.awk $start $end $true_seq_name-log.$ordname-matrix > $true_seq_name-log-info.$ordname-matrix`;


    }
####CREATE DATA STRUCTURE CONTAINING MATRIX OF INTEREST

 	open (PROF,"<$true_seq_name-log-info.$ordname-matrix");
 	while (<PROF>) {
 	    next if m/^#/;
 	    last if m/^\s/;
 	    last if m/^[^\d]/;
 	    chomp;
 	    my @e = split;
 	    push @prof, \@e;
 	}
 	close PROF;
       
     
 	$prof_len = $end - $start + 1;

       chomp ($pwmstring);
#print STDERR "test: $startusrdon or $startusracc or $startusrsta or $startusrbra: 0 or 1 or 0 or 0\n";        
if ($userprofile ) { 
        print SOUT "\nuser-selected $pwmstring profile: start: $start end: $end / $prof_len (profile length)\n\n";
 	print STDERR "\nuser-selected $pwmstring profile: start: $start end: $end / $prof_len (profile length)\n\n";

} else {

        print SOUT "\nautomatically selected $pwmstring profile: start: $start end: $end / $prof_len (profile length)\n\n";
 	print STDERR "\nautomatically selected $pwmstring profile: start: $start end: $end / $prof_len (profile length)\n\n";

}


        return (\@prof,$prof_len,$offset,$start,$end);

#unlink $tempinfolog;
#unlink "$true_seq_name-log-info.$ordname-matrix";
#unlink "$true_seq_name-log-info-pre.$ordname-matrix";

     } ####END GETKMATRIX FUNCTION (Splice site an Start codon PWMs)

     sub numerically { $a <=> $b }

     sub fold4fasta {
	 my $seq = shift;
	 my $foldedseq = "";
	
	 for (my $i=0;$i<length($seq);$i+=60) {
	     my $s = substr($seq,$i,60);

	     $foldedseq = $foldedseq . $s . "\n";
	     #print STDERR $foldedseq;
	 }

	 return $foldedseq;
     }



#Optimize parameter file
 sub OptimizeParameter {

     my($gpfa,$gpgff,$newparam,$branchswitch,$prof_len_bra,$fxdbraoffset,$branchmatrix,$IeWF,$deWF,$FeWF,$IoWF,$doWF,$FoWF,$iMin,$dMin,$fMin,$iAccCtx,$dAccCtx,$fAccCtx) =  @_;
     my @evaluation_total = ();
     my $IeWFini = $IeWF;
     my $IoWFini = $IoWF;
     my $iMinini = $iMin;
     my $iAccCtxini = $iAccCtx;

     if (!$branchswitch) {
     print STDERR "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\n\n";
     print SOUT "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\n\n";
 }else {
     print STDERR "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\nMinBranch range : $iMin to $fMin\nAccCtx range : $iAccCtx to $fAccCtx\n\n";
     print SOUT "\neWF range : $IeWF to $FeWF\noWF range : $IoWF to $FoWF\nMinBranch range : $iMin to $fMin\nAccCtx range : $iAccCtx to $fAccCtx\n\n";

 }

     for ($IeWF=$IeWFini;$IeWF<=$FeWF;$IeWF += $deWF) {
	 print STDERR "eWF: $IeWF\noWF: ";
    
	 for ($IoWF=$IoWFini;$IoWF<=$FoWF;$IoWF += $doWF) {
	     if (!$branchswitch){
	     print STDERR "$IoWF  ";
	 }
	     if ($branchswitch) {
		 print STDERR "$IoWF\nMinDist:  ";
	       for ($iMin=$iMinini;$iMin<=$fMin;$iMin += $dMin) {
		   print STDERR "$iMin\nAccCtx:  ";

		   for ($iAccCtx=$iAccCtxini;$iAccCtx<=$fAccCtx;$iAccCtx += $dAccCtx) {
		        print STDERR "$iAccCtx  ";
		       my $param = Geneid::Param->new();
		       $param->readParam("$results/$species.geneid.param");

		       for (my $i = 0;$i < $param->numIsocores ; $i++) {
			   if (!defined @{$param->isocores}[$i]->Exon_weights([$IeWF,$IeWF,$IeWF,$IeWF])) {
			       die "error in setting exon weights\n";
			   }
			   if (!defined @{$param->isocores}[$i]->Exon_factor([$IoWF,$IoWF,$IoWF,$IoWF])) {
			  # if (!defined @{$param->isocores}[$i]->Exon_factor([0.33,$bestIoWF,$bestIoWF,0.33])) {
			       die "error in setting exon weights\n";
			   }
			   if (!defined @{$param->isocores}[$i]->Site_factor([1-$IoWF,1-$IoWF,1-$IoWF,1-$IoWF])) {
			  # if (!defined @{$param->isocores}[$i]->Site_factor([0.45,1-$bestIoWF,1-$bestIoWF,0.45])) {
			   die "error in setting exon weights\n";
			   }
			   if (!defined @{$param->isocores}[$i]->set_profile('Branch_point_profile',$prof_len_bra,$fxdbraoffset,-50,0,0,1,$iAccCtx,$iMin,0,0,$branchmatrix)){die "error in setting profile\n";}
		       }


		       $param->writeParam("$results/$species.geneid.param.temp");
            
			print STDERR "geneid -GP $results/${newparam}.temp $gpfa \| egrep -vw 'exon' \|  gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' > $TMP/Predictions.${newparam}.gff";
	     `geneid -GP $results/${newparam}.temp $gpfa | egrep -vw 'exon' | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' > $TMP/Predictions.${newparam}.gff`;
	

	     my  @evaluation_output  = split " ",  `evaluation -sta $TMP/Predictions.${newparam}.gff $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $IeWF, $iAccCtx, $iMin, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}'`; 



 push (@evaluation_total,\@evaluation_output);
	    
		   }#$iAccCtx for loop
		   unless ($iMin == $fMin) {print STDERR "\nMinDist:  "};
		    $iAccCtx = $iAccCtxini;

	       }#$iMin for loop
	       unless ($IoWF == $FoWF) {print STDERR "\noWF:  "};
	       $iAccCtx = $iAccCtxini;
	       $iMin=$iMinini;
	       
	       
	   }#if $branchswitch
	     elsif (!$branchswitch) {


my $param = Geneid::Param->new();
		       $param->readParam("$results/$species.geneid.param");

		       for (my $i = 0;$i < $param->numIsocores ; $i++) {
			   if (!defined @{$param->isocores}[$i]->Exon_weights([$IeWF,$IeWF,$IeWF,$IeWF])) {
			       die "error in setting exon weights\n";
			   } 
			  if (!defined @{$param->isocores}[$i]->Exon_factor([$IoWF,$IoWF,$IoWF,$IoWF])) {
		 #   if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$IoWF,$IoWF,0.4])) {
			       die "error in setting exon weights\n";
			   } 
			   if (!defined @{$param->isocores}[$i]->Site_factor([1-$IoWF,1-$IoWF,1-$IoWF,1-$IoWF])) {
		#	  if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$IoWF,1-$IoWF,0.55])) {
				die "error in setting exon weights\n";
			   } 
		       
		       
		       }

		       $param->writeParam("$results/$species.geneid.param.temp");


	     `geneid -GP $results/${newparam}.temp $gpfa | egrep -vw 'exon' | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' > $TMP/Predictions.${newparam}.gff`;
	

	     my  @evaluation_output  = split " ",  `evaluation -sta $TMP/Predictions.${newparam}.gff $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $IoWF, $IeWF, \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}'`; 



 push (@evaluation_total,\@evaluation_output);





	     } #elseif
	       

	 }

	# $iAccCtx = $iAccCtxini;
	# $iMin=$iMinini;
	 $IoWF = $IoWFini;
	 print STDERR "\n";

     }



 return \@evaluation_total;




 } # end sub optimize parameter file


sub BuildOptimizedParameterFile{

 my($evalarray,$branchswitch,$prof_len_bra,$fxdbraoffset,$branchmatrix) =  @_;
 my @sortedeval = ();
 my @evaluationinit= ();
 my $bestIoWF = "";
 my $bestIeWF = "";
 my $bestMin = "";
 my $bestAcc = "";

 if (!$branchswitch){

 @sortedeval = sort sorteval @$evalarray;

 $bestIoWF = $sortedeval[0][0]; #0.2
 $bestIeWF = $sortedeval[0][1]; #-3.5

print STDERR "\nBest performance obtained using IoWF: ".$sortedeval[0][0]." and IeWF: ".$sortedeval[0][1]."\n";
print SOUT "\nBest parameter file performance obtained using oWF: ".$sortedeval[0][0]." and eWF: ".$sortedeval[0][1]."\n";

#INITIALIZE ARRAY WITH EVALUATION PARAMETERS
@evaluationinit= (qw(oWF eWF SN SP CC SNe SPe SNSP SNg SPg SNSPg raME raWE));

print STDERR "Sorted performance results (Three best performance estimates) for different values of oWF and eWF:\n".join("\t",@evaluationinit),"\n";
print SOUT "\nSorted performance results (best to worst) for different values of oWF and eWF:\n\n".join("\t",@evaluationinit),"\n";

foreach my $eval_ref (@sortedeval){

print SOUT join("\t",@$eval_ref),"\n";

}

###FOUR BEST PERFORMANCE RESULTS AFTER OPTIMIZATION
for (my $i=0; $i<=2;$i++)
  {
print STDERR join("\t",@{$sortedeval[$i]}),"\n";
  }
############

#BUILD BEST PERFORMANCE (OPTIMIZED) PARAMETER FILE (USE VALUES OBTAINED ABOVE)

$param = Geneid::Param->new();
	     $param->readParam("$results/$species.geneid.param");



	     for (my $i = 0;$i < $param->numIsocores ; $i++) {
		 if (!defined @{$param->isocores}[$i]->Exon_weights([$bestIeWF,$bestIeWF,$bestIeWF,$bestIeWF])) {
		     die "error in setting exon weights\n";
		 } 
		 if (!defined @{$param->isocores}[$i]->Exon_factor([$bestIoWF,$bestIoWF,$bestIoWF,$bestIoWF])) {
	#	 if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$bestIoWF,$bestIoWF,0.4])) {
		 die "error in setting exon weights\n";
		 } 
		 if (!defined @{$param->isocores}[$i]->Site_factor([1-$bestIoWF,1-$bestIoWF,1-$bestIoWF,1-$bestIoWF])) {
	#	 if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$bestIoWF,1-$bestIoWF,0.55])) {
		 die "error in setting exon weights\n";
		 } 
	     }
#write new parameter file (optimized)
$param->writeParam("$results/$species.geneid.optimized.param");

print STDERR "\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";
print SOUT "\n\nNew optimized parameter file named: $species.geneid.optimized.param \n";

return [$bestIeWF,$bestIoWF,0,0,\@evaluationinit];

} elsif ($branchswitch){

 @sortedeval = sort sortevalbranch @$evalarray;

 $bestIoWF = $sortedeval[0][0]; #0.2
 $bestIeWF = $sortedeval[0][1]; #-3.5
 $bestAcc = $sortedeval[0][2];
 $bestMin = $sortedeval[0][3];
 


print STDERR "\nBest performance obtained using IoWF: ".$sortedeval[0][0]." , IeWF: ".$sortedeval[0][1].", MinimalBranchDist: ".$bestMin.", Optimal Branch Context: ".$bestAcc."\n";
print SOUT "best performance obtained using IoWF: ".$sortedeval[0][0]." , IeWF: ".$sortedeval[0][1].", MinimalBranchDist: ".$bestMin.", Optimal Branch Context: ".$bestAcc."\n";

#INITIALIZE ARRAY WITH EVALUATION PARAMETERS
@evaluationinit= (qw(oWF eWF AccCtx MinD SN SP CC SNe SPe SNSP SNg SPg SNSPg raME raWE));

print STDERR "(Sorted performance results (Three best performance estimates) for different values of oWF, eWF, AccCtx and MinD)\n\n".join("\t",@evaluationinit),"\n";
print SOUT "(Sorted performance results (best to worst) for different values of oWF, eWF, AccCtx and MinD)\n\n".join("\t",@evaluationinit),"\n";

foreach my $eval_ref (@sortedeval){

print SOUT join("\t",@$eval_ref),"\n";

}

###THREE BEST PERFORMANCE RESULTS AFTER OPTIMIZATION
for (my $i=0; $i<=2;$i++)
  {
print STDERR join("\t",@{$sortedeval[$i]}),"\n";
  }
############

#BUILD BEST PERFORMANCE (OPTIMIZED) PARAMETER FILE (USE VALUES OBTAINED ABOVE)

$param = Geneid::Param->new();
	     $param->readParam("$results/$species.geneid.param");



	     for (my $i = 0;$i < $param->numIsocores ; $i++) {
		 if (!defined @{$param->isocores}[$i]->Exon_weights([$bestIeWF,$bestIeWF,$bestIeWF,$bestIeWF])) {
		     die "error in setting exon weights\n";
		 } 
	         if (!defined @{$param->isocores}[$i]->Exon_factor([$bestIoWF,$bestIoWF,$bestIoWF,$bestIoWF])) {
	       #   if (!defined @{$param->isocores}[$i]->Exon_factor([0.4,$bestIoWF,$bestIoWF,0.4])) {
		     die "error in setting exon weights\n";
		 } 
		 if (!defined @{$param->isocores}[$i]->Site_factor([1-$bestIoWF,1-$bestIoWF,1-$bestIoWF,1-$bestIoWF])) {
		#  if (!defined @{$param->isocores}[$i]->Site_factor([0.55,1-$bestIoWF,1-$bestIoWF,0.55])) {
		    die "error in setting exon weights\n";
		 } 
		 if (!defined @{$param->isocores}[$i]->set_profile('Branch_point_profile',$prof_len_bra,$fxdbraoffset,-50,0,0,1,$bestAcc,$bestMin,0,0,$branchmatrix)){die "error in setting profile\n";}
	     }
#write new parameter file (optimized)
$param->writeParam("${results}/$species.geneid.optimized.param");

print STDERR "\nNew optimized parameter file named: $species.geneid.optimized.param \n";
print SOUT "\nNew optimized parameter file named: $species.geneid.optimized.param \n";


return [$bestIeWF,$bestIoWF,$bestAcc,$bestMin,\@evaluationinit];

} #if branch switch



}



sub EvaluateParameter {

     my($gpfa,$gpgff,$newparam) =  @_;
    
     
	     `geneid -GP $results/$newparam $gpfa | egrep -vw 'exon' | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' > $results/test.predictions.${newparam}.gff`;
	

	     my  @evaluation_test  = split " ",  `evaluation -sta $results/test.predictions.${newparam}.gff $gpgff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", \$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}'`; 


     return \@evaluation_test;




 } # evaluate parameter function

###JACKNIFE SUB

sub runJacknife {

    my($list4jkf,$seqstofill,$grpsiz,$outacceptortbl,$outdonortbl,$outstarttbl,$bckgrnd,$outcds,$outintron,$bestIoWF,$bestIeWF,$shortintron,$longintron,$minintergenic,$maxintergenic,$gptraintbl,$gptraingff,$startdonor,$enddonor,$startacceptor,$endacceptor,$startstart,$endstart,$branchsw,$fullengthbranchtbl,$startbranch,$endbranch,$bestAcc,$bestMin) =  @_;
    my $grepBranches = "";
    my $branchexcept = "";
    my $count = 1;
    my $geneidout = "";
    #my $cycles = 402;
    my $cycles = int (scalar(@$list4jkf) / $grpsiz); #610/10
    my $status = "";
  #  my $tempgroup4jckf = "";
    my $actualannots = scalar(@$list4jkf) - $seqstofill;

    print STDERR "For the 10-fold cross-validation method we will remove groups of $grpsiz sequences from the ".$actualannots." training set annotations, train on the remaining annotations and run geneid evaluate the performance on the left-out sequences (performed $cycles times)\n\n";
    print SOUT "For the 10-fold cross-validation method we will remove groups of $grpsiz sequences from the ".$actualannots." training set annotations, train on the remaining annotations and run geneid evaluate the performance on the left-out sequences (performed $cycles times)\n\n";

###############################
#for loop to go through seqs
#########################3    
    for (my $i = 0; $i < scalar(@$list4jkf); $i += $grpsiz) {
  
	print STDERR "Extracting the set $count of sequences for 10x cross validation training and evaluation ($count out of $cycles)\n";

	my @seqsjk = @$list4jkf[$i .. $i+$grpsiz-1];
	my $listseqsjacknife = join ("\n",@seqsjk);
#        print STDERR "group1:".$listseqsjacknife."\n\n";
	# my $listseqsjacknife4eval =  "$listseqsjacknife"."#$";
	
       my $group4jckf = $species."_seqs_leftout";   
	open FOUT, ">$group4jckf";
	print FOUT "$listseqsjacknife";
	close FOUT;
	
	##grepSeqsToEval

	my $seqsaside1 = "";
	open LOCID, "sed -e '/^\$/d' -e 's/\\(.*\\)/\\1\$/g' $group4jckf |";
	while (<LOCID>) {
	    $seqsaside1 .= $_;
	}
	close LOCID;

	my $tempgroup4jckf1 = $species."_seqs_leftout_crossvalidation_1";   
	open FOUT, ">$tempgroup4jckf1";
	print FOUT "$seqsaside1";
	close FOUT;
#		print STDERR " $tempgroup4jckf1\n";
	#Acceptors,Donors,Starts, Noncoding,branches
	
	my $seqsaside2 = "";
	open LOCID, "sed -e '/^\$/d' -e 's/\\(.*\\)/\\1\./g' $group4jckf |";
	while (<LOCID>) {
	    $seqsaside2 .= $_;
	}
	close LOCID;

	my $tempgroup4jckf2 = $species."_seqs_leftout_crossvalidation_2";   
	open FOUT, ">$tempgroup4jckf2";
	print FOUT "$seqsaside2";
	close FOUT;
#	         print STDERR " $tempgroup4jckf2\n";

	#Coding
	
	my $seqsaside3 = "";
	open LOCID, "sed -e '/^\$/d' -e 's/\\(.*\\)/\\1\_/g' $group4jckf |";
	while (<LOCID>) {
	    $seqsaside3 .= $_;
	}
	close LOCID;

	my $tempgroup4jckf3 = $species."_seqs_leftout_crossvalidation_3";   
	open FOUT, ">$tempgroup4jckf3";
	print FOUT "$seqsaside3";
	close FOUT;
#	                 print STDERR " $tempgroup4jckf3\n";
	

	#CREATE BLANK TEMP PARAMETER FILE FOR JACKNIFE
	my $paramtemp = Geneid::Param->new($species);
	#set isochores to 1
	$paramtemp->numIsocores(1);
	$paramtemp->isocores([Geneid::Isocore->new()]);
	
#	print SDTERR $outacceptortbl."\n";
	#####PRE EXTRACT SEQUENCES FROM TRAINING SET
	my $grepAcceptor="egrep -vf $tempgroup4jckf2 $outacceptortbl  > $TMP/tmp_Acceptorsminus";

	my $grepDonor="egrep -vf $tempgroup4jckf2 $outdonortbl  > $TMP/tmp_Donorsminus";
	
	my $grepStarts="egrep -vf $tempgroup4jckf2 $outstarttbl  > $TMP/tmp_Startsminus";
	
	my $grepMarkovCoding="egrep -vf $tempgroup4jckf3 $outcds  > $TMP/tmp_Codingminus";
	
	my $grepMarkovNonCoding="egrep -vf $tempgroup4jckf2 $outintron  > $TMP/tmp_NonCodingminus";
	
#	print STDERR "$gptraintbl $tempgroup4jckf1";
	
	my $grepSeqsToEval="gawk '{print \$2,\$1}' $gptraintbl | egrep -f $tempgroup4jckf1 - | gawk '{print \$2,\$1}' - > $TMP/tmp_SeqsToEval";
		
	if ($branchsw){
        $grepBranches="egrep -vf $tempgroup4jckf2 $fullengthbranchtbl  > $TMP/tmp_Branchesminus";


    }
	$status= system $grepAcceptor;

	if ($status != 0) {
	  print STDERR "FATAL ERROR !!! Unsuccessful command :\n  $grepAcceptor" && die;
	};
#	print STDERR "$grepAcceptor\n";
 
	$status= system $grepDonor;

	if ($status != 0) {
	  print STDERR "FATAL ERROR !!! Unsuccessful command :\n $grepDonor" && die;
	};
#	print STDERR "$grepDonor\n";

  
	$status= system $grepStarts;

	if ($status != 0) {
	  print STDERR "FATAL ERROR !!! Unsuccessful command :\n  $grepStarts" && die;
	};
#	print STDERR "$grepStarts\n";

if ($branchsw){	
	$status= system $grepBranches;

	if ($status != 0) {
	  print STDERR "FATAL ERROR !!! Unsuccessful command :\n  $grepBranches" && die;
	};

    }
	$status= system $grepMarkovCoding;

	if ($status != 0) {
	  print STDERR "FATAL ERROR !!! Unsuccessful command :\n $grepMarkovCoding" && die;
	};
#	print STDERR "$status $grepMarkovCoding\n";
  
	$status= system $grepMarkovNonCoding;

	if ($status != 0) {
	   print STDERR "FATAL ERROR !!! Unsuccessful command :\n  $grepMarkovNonCoding" && die;
	};
#	print STDERR "$grepMarkovNonCoding\n";

	$status= system $grepSeqsToEval;

	if ($status != 0) {
	   print STDERR "FATAL ERROR !!! Unsuccessful command :\n $grepSeqsToEval" && die;
	};
#	print SDTERR "$status $grepSeqsToEval";

###store taken out sequences in the following variables:

	my $acceptorexcept = "$TMP/tmp_Acceptorsminus";
	my $donorexcept = "$TMP/tmp_Donorsminus";
	my $startexcept = "$TMP/tmp_Startsminus";
	my $codingexcept = "$TMP/tmp_Codingminus";
	my $noncodingexcept = "$TMP/tmp_NonCodingminus";
	my $seqstoevaltbl = "$TMP/tmp_SeqsToEval";
if ($branchsw){
        $branchexcept = "$TMP/tmp_Branchesminus";
    
    }

	###CONVERT TO FASTA (EVALUATE SET)
	my $tempgp_fa_eval_jkf = $species.".gp.eval.crossvalidation.fa";   
        TblToFasta($seqstoevaltbl,$tempgp_fa_eval_jkf);

	



	#####EXTRACT SEQUENCES FROM JACKNIFE TRAINING SET

	########
	#########
	#get donor site statistics
	#########
	#########
	my $order = "0";
	my $numbersites = `wc $donorexcept | gawk '{print \$1}'`;
	chomp $numbersites;
	my $donoffset = "30"; #position before intron (last of exon (31) -1 for offset)

	if ($numbersites > 1200) {

	    $order = "1";

	} elsif ($numbersites <= 1200) {

	    $order = "0";
	}

	print STDERR "$startdonor/$enddonor/$prof_len_don There are $numbersites donor sites, enough for a matrix of order $order (Jacknife) set $count \n";

	my ($donormatrix,$prof_len_don,$fxddonoffset,$s,$e) = getKmatrix($donorexcept,$bckgrnd,$order,$donoffset,1,0,0,0,$startdonor,$enddonor,1);
	if (!defined @{$paramtemp->isocores}[0]->set_profile('Donor_profile',$prof_len_don,$fxddonoffset,$cutoff,$order,0,1,0,0,0,0,$donormatrix)) {
	    die "error in setting profile\n";
	}


	########($true_seqs,$false_seqs,$order,$offset,$donor,$accept,$star,$branch,$start,$end,$jacknife)
	#########
	#get acceptor site statistics
	#########
	#########
	$order = "0";
	$numbersites = `wc $acceptorexcept | gawk '{print \$1}'`;
	chomp $numbersites;
	my $accoffset = "30"; #position after intron (first of exon (31) -1 for offset)

	if ($numbersites > 1200) {

	    $order = "1";

	} elsif ($numbersites <= 1200) {

	    $order = "0";
	}

	print STDERR "$startacceptor/$endacceptor/$prof_len_acc:There are $numbersites acceptor sites, enough for a matrix of order $order (Jacknife) set $count \n";

	my ($acceptormatrix,$prof_len_acc,$fxdaccoffset,$st,$en) = getKmatrix($acceptorexcept,$bckgrnd,$order,$accoffset,0,1,0,0,$startacceptor,$endacceptor,1);
	if (!defined @{$paramtemp->isocores}[0]->set_profile('Acceptor_profile',$prof_len_acc,$fxdaccoffset,$cutoff,$order,0,1,0,0,0,0,$acceptormatrix)) {
	    die "error in setting profile\n";
	}

	#########my ($true_seqs,$false_seqs,$order,$offset,$donor,$accept,$star,$branch,$start,$end,$jacknife)
	#########
	#get start site statistics
	#########
	#########
	$order = "0";
	$numbersites = `wc $startexcept | gawk '{print \$1}'`;
	chomp $numbersites;
	my $staoffset = "30"; #before first position of the exon (31)minus 1 for offset)

	if ($numbersites > 5500) {

	    $order = "2";

	} elsif ($numbersites <= 5500) {

	    $order = "0";
	}

	print STDERR "There are $numbersites start sites, enough for a matrix of order $order (Jacknife) set $count \n";

	my ($startmatrix,$prof_len_sta,$fxdstaoffset,$sta,$enj) = getKmatrix($startexcept,$bckgrnd,$order,$staoffset,0,0,1,0,$startstart,$endstart,1);

	####write to parameter file
	if (!defined @{$paramtemp->isocores}[0]->set_profile('Start_profile',$prof_len_sta,$fxdstaoffset,$cutoff,$order,0,1,0,0,0,0,$startmatrix)) {
	    die "error in setting profile\n";
	}
	#############################my ($true_seqs,$false_seqs,$order,$offset,$donor,$accept,$star,$branch,$start,$end,$jacknife) = @_;


	if ($branchsw) {

	#########
	#get branch site statistics
	#########
	#########
	$order = "0";
	$numbersites = `wc $branchexcept | gawk '{print \$1}'`;
	chomp $numbersites;
	my $braoffset = "32"; #before first position of the exon (35)minus 1 for offset)

	print STDERR "There are $numbersites branch sites, enough for a matrix of order $order (Jacknife) set $count \n";

	my ($branchmatrix,$prof_len_bra,$fxdbraoffset,$startbranch,$endbranch) = getKmatrix($fullengthbranchtbl,$bckgrnd,$order,$braoffset,0,0,0,1,$startbranch,$endbranch,1);

	####write to parameter file
	if (!defined @{$paramtemp->isocores}[0]->set_profile('Branch_point_profile',$prof_len_bra,$fxdbraoffset,"-50",$order,0,1,$bestAcc,$bestMin,0,0,$branchmatrix)) {
	    die "error in setting profile\n";
	}
	#############################

    }


	###DERIVE INITIAL/TRANSITION MARKOV MODEL JACKNIFE


	#print STDERR "\nDeriving markov model (Jacknife method): set $count \n";

	my ($markovini,$markovtrans,$totalcoding,$totalnoncoding,$markovmodel) = @{deriveCodingPotential($codingexcept,$noncodingexcept)};

	#add markov matrices to the parameter file 
	if (!defined @{$paramtemp->isocores}[0]->Markov_order($markovmodel)) {
	    die "error in setting Markov_order\n";
	}
	if (!defined @{$paramtemp->isocores}[0]->Markov_Initial_probability_matrix($markovini)) {
	    die "error in setting Markov_Initial_probability_matrix\n";
	}
	if (!defined @{$paramtemp->isocores}[0]->Markov_Transition_probability_matrix($markovtrans)) {
	    die "error in setting Markov_Transition_probability_matrix\n";
	}

	###WRITE PRELIMINARY PARAMETER FILE FOR JACKNIFE
	#$paramtemp->writeParam("$species.geneid.param.tmp");

	#my $newparamtemp = "$species.geneid.param.tmp";

	#$paramtemp = Geneid::Param->new();
	#	     $param->readParam("$newparamtemp");

	for (my $i = 0;$i < $paramtemp->numIsocores ; $i++) {
	    if (!defined @{$paramtemp->isocores}[$i]->Exon_weights([$bestIeWF,$bestIeWF,$bestIeWF,$bestIeWF])) {
		die "error in setting exon weights\n";
	    } 
	    if (!defined @{$paramtemp->isocores}[$i]->Exon_factor([$bestIoWF,$bestIoWF,$bestIoWF,$bestIoWF])) {
	  #   if (!defined @{$paramtemp->isocores}[$i]->Exon_factor([0.4,$bestIoWF,$bestIoWF,0.4])) {
		die "error in setting exon weights\n";
	    } 
	   if (!defined @{$paramtemp->isocores}[$i]->Site_factor([1-$bestIoWF,1-$bestIoWF,1-$bestIoWF,1-$bestIoWF])) {
	  #  if (!defined @{$paramtemp->isocores}[$i]->Site_factor([0.55,1-$bestIoWF,1-$bestIoWF,0.55])) {
		die "error in setting exon weights\n";
	    } 
	} # for number of isochores

	#print STDERR "bestIoWF jacknife: $bestIoWF\n";
        #Open gene model object
	$paramtemp->geneModel(Geneid::GeneModel->new());
	$paramtemp->geneModel->useDefault;
	$paramtemp->geneModel->intronRange($shortintron,$longintron);
	#$paramtemp->geneModel->intergenicRange(200,'Infinity');
	$paramtemp->geneModel->intergenicRange($minintergenic,$maxintergenic);
	###WRITE PRELIMINARY PARAMETER FILE FOR JACKNIFE
	$paramtemp->writeParam("$species.geneid.jacknife.param.temp.$count");
	my $newparamtemp = "$species.geneid.jacknife.param.temp.$count";

	####RUN GENEID ON THE LEFT OUT SEQUENCES:
	print STDERR "\nRunning geneid on set number $count\n\n\n";

	open LOCID, "geneid -GP $results/$newparamtemp $tempgp_fa_eval_jkf | egrep -wv 'exon' | gawk 'NR>5 {if (\$2==\"Sequence\") print \"\#\$\"; if (substr(\$1,1,1)!=\"\#\") print }' |";
	while (<LOCID>) {
	    $geneidout .= $_;
	}
	close LOCID;


	$geneidout .= '#$'."\n";

$count++;
    
    } #end of for loop				###FOR LOOP: GO THROUGH GROUP OF SEQS


   ###### ###EVALUATE GENEID AGAINST ANNOTATIONS

    $temp_jkf_geneid = $species.".geneid_crossvalidation";   
    open FOUT, ">$temp_jkf_geneid";
    print FOUT "$geneidout";
    close FOUT;

if ($branchsw) {

    my  @evaluation_jacknife  = split " ",  `evaluation -sta $temp_jkf_geneid $gptraingff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\%6.2f \%6.2f \%6.2f \%6.2f\\n\", $bestIoWF,$bestIeWF,$bestAcc,$bestMin,\$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}'`; 

    return \@evaluation_jacknife;

} elsif (!$branchsw){
    
    my  @evaluation_jacknife  = split " ",  `evaluation -sta $temp_jkf_geneid $gptraingff | tail -2 | head -1 |  gawk '{printf \"\%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f \%6.2f\\n\", $bestIoWF,$bestIeWF,\$1, \$2, \$3, \$4, \$5, \$6, \$9, \$10, \$11, \$7, \$8}'`; 

    return \@evaluation_jacknife;

}


#unlink $temp_jkf_geneid;

   }				#SUB JACKNIFE

sub WriteStatsFile {


my($species,$outintron,$outcds,$outgff,$inframe,$inframeeval,$seqsused,$totalnoncandon,$totalnoncanacc,$totalnoncansta,$markovmodel,$totalcoding,$totalnoncoding,$stdo,$endo,$stac,$enac,$stst,$enst,$stbr,$enbr,$branchsw,$useallseqs) =  @_;

###OBTAIN GENE MODEL SET STATISTICS

#Open gene model object

$param->geneModel(Geneid::GeneModel->new());
$param->geneModel->useDefault;
###########
#open SOUT,">$sout";

my $avgintron = "";
my $sdintron = "";

my @intronlength = `gawk '{print length(\$2)}' $outintron`;
my ($mean, $st) = average(\@intronlength); 
#print STDERR "INTRON: $mean, $st\n";

###INTRON MIN AND MAX LENGTH

#if (!$useextdata) {

my @intronlist = `gawk '{print length(\$2)}' $outintron | sort -n`;

my $totintrons = scalar(@intronlist);
my @intronlen =();


my $intr = "";
for (my $i=0;$i<=scalar(@intronlist)-1;$i++) {

   $intr = $intronlist[$i];
   chomp $intr;
   push (@intronlen, $intr);
}

my @slice1 = @intronlen[0..5];
my @slice2 = @intronlen[(scalar(@intronlen) - 5)..(scalar(@intronlen) - 1)];

my $shortintron = int($intronlist[0] - $intronlist[0]*(0.25));
chomp $shortintron;
#my $longintron =  $intronlist[$totintrons - 1];
my $longintron = int($mean + ($st * 2) > 6000 ? 25000 : $mean + ($st * 2));
chomp $longintron;

###} #execute if not using external data 


#my $answer = "";
#do {
 #            print "\nDo you want to use this automatically selected range of introns ($shortintron to $longintron) in the geneid gene model (yes/no)?\n(Note that the 5 smallest introns were found to be: ".join(", ",@slice1)." nucleotides long and the 5 longest introns: ".join(", ",@slice2)." bases in length)\n";
  #	    $answer=<STDIN>;
#	     chomp $answer;
 #	} while ($answer !~ /^(yes|y)|(n|no)$/i);

if (!$useextdata ||  $shortintronusr==0) {

print "\n The selected range of introns in the geneid gene model was automatically set to be ($shortintron to $longintron) \n(Note that the 5 smallest introns were found to be: ".join(", ",@slice1)." nucleotides long and the 5 longest introns: ".join(", ",@slice2)." bases in length)\n";

} #execute if not using external data 
elsif ($useextdata && $shortintronusr>0) {

    #my $shortintronaut=$shortintron;
    #my $longintronaut=$longintron;


print "\n The user-selected range of introns in the geneid gene model was set to be ($shortintronusr to $longintronusr) rather than the automatically selected $shortintron to $longintron \n(Note that the 5 smallest introns were found to be: ".join(", ",@slice1)." nucleotides long and the 5 longest introns: ".join(", ",@slice2)." bases in length)\n";

 }


#print SOUT  "\n The selected range of introns in the geneid gene model was set to be ($shortintron to $longintron) \n(Note that the 5 smallest introns were found to be: ".join(", ",@slice1)." nucleotides long and the 5 longest introns: ".join(", ",@slice2)." bases in length)\n";

=head
if ($answer =~ /^(no|n)$/i) {
 my $sline = "";
 my $eline = "";
 	    do {
		print STDERR "\nType new intron minimum length boundary: ";
		$sline = readline(STDIN);
 	    } while ($sline !~/(\d+)/);
 	    $shortintron = $1;
 	    do {
		print STDERR "\nType maximum intron length boundary: ";
		$eline = readline(STDIN);
 	    } while ($eline !~/(\d+)/ || $eline <= $sline);
 	    $longintron = $1;
}
=cut

my $minintergenic = 200;
my $maxintergenic = 'Infinity';

#do {
 #            print "\nDo you want to use this automatically selected intergenic distance range ($minintergenic to $maxintergenic) in the geneid gene model (yes/no)?\n";
  #	    $answer=<STDIN>;
#	     chomp $answer;
 #	} while ($answer !~ /^(yes|y)|(n|no)$/i);

if (!$useextdata || $minintergenicusr==0) {

print "\nThe automatically selected intergenic distance range was set to be ($minintergenic to $maxintergenic) in the geneid gene model\n";

} elsif ($useextdata && $minintergenicusr>0) {

   #my $minintergenicaut=$minintergenic;
    #my $maxintergenicaut=$maxintergenic;

print "\nThe user-selected intergenic distance range was set to be ($minintergenicusr to $maxintergenicusr) rather than the default ($minintergenic to $maxintergenic) in the geneid gene model\n";


}
 

#print SOUT  "\nThe automatically selected intergenic distance range was set to be ($minintergenic to $maxintergenic) in the geneid gene model\n";

=head
if ($answer =~ /^(no|n)$/i) {
    my $sline = "";
    my $eline = "";
 	    do {
		print STDERR "\nType new minimum intergenic distance: ";
		$sline = readline(STDIN);
 	    } while ($sline !~/(\d+)/);
 	    $minintergenic = $1;
 	    do {
		print STDERR "\nType maximum intergenic distance (type a number or 'Infinity'): ";
		$eline = readline(STDIN);
 	    } while ($eline !~/(\d+|Infinity)/ || $eline <= $sline );
 	    $maxintergenic = $1;

}
=cut

##use shortest and longest intron lengths in gene model of parameter file 
if (!$useextdata || ($minintergenicusr==0 && $shortintronusr==0)) {

$param->geneModel->intronRange($shortintron,$longintron);
$param->geneModel->intergenicRange($minintergenic,$maxintergenic);

} elsif ($useextdata && $minintergenicusr>0 && $shortintronusr>0) {

$param->geneModel->intronRange($shortintronusr,$longintronusr);
$param->geneModel->intergenicRange($minintergenicusr,$maxintergenicusr);

} elsif ($useextdata && $minintergenicusr==0  && $shortintronusr>0) {

$param->geneModel->intronRange($shortintronusr,$longintronusr);
$param->geneModel->intergenicRange($minintergenic,$maxintergenic);

} elsif ($useextdata && $minintergenicusr>0  && $shortintronusr==0) {

$param->geneModel->intronRange($shortintron,$longintron);
$param->geneModel->intergenicRange($minintergenicusr,$maxintergenicusr);

}

###############################

my @CDSGCcontent = `gawk '{print gsub(/[GC]/,".",\$2)/length(\$2)}' $outcds`;
#print STDERR "@CDSGCcontent\n";
my ($meangc,$stgc) = average(\@CDSGCcontent);
#print STDERR "CDS: $meangc $stgc $outintron\n";

my @intronGCcontent = `gawk '{print gsub(/[GC]/,".",\$2)/length(\$2)}' $outintron`;

#print STDERR "@intronGCcontent\n";
my ($meangci,$stgci) = average(\@intronGCcontent);
#print STDERR "intron: $meangci $stgci\n";

my $totexons = `gawk '{print \$9}' $outgff | wc | gawk '{print \$1}' `;
chomp $totexons;

my @exonspergene = `gawk '{print \$9}' $outgff | sort | uniq -c | gawk '{print \$1}'`;

my ($avgex,$stex) = average(\@exonspergene);

my @exonlength = `egrep -v 'Single' $outgff | gawk '{len=\$5-\$4;print len}' - | sort`;

my ($avgle,$stle) = average(\@exonlength);

my $singlegenes = `egrep -c '(Single)' $outgff`; 
chomp $singlegenes;

#print SOUT "GENE MODEL STATISTICS FOR $species\n\n";

print SOUT "\n\nA subset of $totalseqs4training sequences (randomly chosen from the $total_seqs gene models) was used for training\n\n";

if (!$useallseqs) {
print SOUT "\n$seqsused gene models (80 % of total) were used for training and $gffseqseval annotations (20 % of total) set aside for evaluation (randomly)\n\n";
} else {print SOUT "\n$total_seqs gene models were used for both training and evaluation\n\n"; }

if (!$useallseqs){
    print SOUT "\n$inframe of the gene models translate into proteins with in-frame stops within the training set and $inframeeval in the evaluation set (seqs removed).\n\n"; } else { print SOUT "$inframe of the gene models translate into proteins with in-frame stops within the training set.\n\n"; }

print SOUT "There are $totalnoncandon non-canonical donors as part of the training set\n\n";
print SOUT "There are $totalnoncanacc non-canonical acceptors as part of the training set\n\n";
print SOUT "There are $totalnoncansta non-canonical start sites as part of the training set\n\n";
print SOUT "These gene models correspond to $totalcoding coding bases and $totalnoncoding non-coding bases\n\n";
print SOUT "Deriving a markov model for the coding potential of order $markovmodel\n\n";
print SOUT "The intronic sequences extracted from the gene models have an average length of $mean, with $st of SD\n\n";

if (!$useextdata || $shortintronusr<="0") {
print SOUT "Geneid can predict gene models having introns with a minimum length of $shortintron nucleotides and a maximum of $longintron bases (boundaries used in gene model were automatically selected) \n\n";
} elsif ($useextdata && $shortintronusr>"0") {
print SOUT "Geneid can predict gene models having introns with a minimum length of $shortintronusr nucleotides and a maximum of $longintronusr bases (boundaries used in gene model which were user-selected) \n\n";
}


if (!$useextdata || $minintergenicusr<="0") {
print SOUT "The minimum (automatically selected) intergenic distance was set to $minintergenic nucleotides whereas the maximum was set to $maxintergenic (boundaries used in gene model) \n\n";
} elsif ($useextdata && $minintergenicusr>"0") {
    print SOUT "The minimum (user selected) intergenic distance was set to $minintergenicusr nucleotides whereas the maximum was set to $maxintergenicusr (boundaries used in gene model) \n\n"; }

print SOUT "The GC content of the exonic and intronic sequences is $meangc (SD $stgc) and $meangci (SD $stgci) respectively \n\n";
print SOUT "The gene models used for training contain $totexons exons \n\n";
print SOUT "The gene models average $avgex exons per gene (SD $stex)\n\n";
print SOUT "The average length of the exons (non-single) in the training set gene models is $avgle (SD $stle)\n\n";

if (!$useallseqs) {
print SOUT "The training set includes $singlegenes single-exon genes (out of $seqsused ) gene models\n\n";
} else {print SOUT "The training set includes $singlegenes single-exon genes (out of $total_seqs) gene models\n\n"; }


print SOUT "The donor site profile chosen by the user spans ".($endo-$stdo+1)." nucleotides: position $stdo to $endo\n";
print SOUT "The acceptor site profile chosen by the user spans ".($enac-$stac+1)." nucleotides: position $stac to $enac\n";
print SOUT "The start site profile chosen by the user spans ".($enst-$stst+1)." nucleotides: position $stst to $enst\n";
	if ($branchsw) {
print SOUT "The branch site profile chosen by the user spans ".($enbr-$stbr+1)." nucleotides: position $stbr to $enbr\n";
	}

if (!$useextdata || ($minintergenicusr=="0" && $shortintronusr=="0")) {

return ($shortintron,$longintron,$minintergenic,$maxintergenic);

} 

elsif ($useextdata && $shortintronusr=="0") {

return ($shortintron,$longintron,$minintergenicusr,$maxintergenicusr);

       }   

elsif ($useextdata && $minintergenicusr=="0") {

return ($shortintronusr,$longintronusr,$minintergenic,$maxintergenic);

       }
	 
}


sub average {
    my ($sequences) = @_;
    my $sum = 0;
    my $total = 0;
    my ($mean, $st);

    foreach my $seq (@$sequences) {
	$sum += $seq;
	$total++;
    }
    
    $mean = $sum / $total;
    $mean = sprintf("%.3f",$mean);

    $sum = 0;

    foreach my $seq (@$sequences) {
	$sum += ($seq - $mean) * ($seq - $mean);

    }

    $st = sqrt($sum / $total);
    $st = sprintf("%.3f",$st);

    return($mean,$st); 
    
}



sub predictPlotgff2ps {

    my ($paramopt,$gpfa,$gpgff,$gplen,$tempjkf_geneid) = @_;


    print STDERR "\nRunning geneid on fastas\n\n";
 
    my $geneidall = "";
    my $gff2psplots =  $plotsdir."$species.gff2ps.prediction.plots.ps";
  

     open LOCID, "geneid -GP $results/$paramopt $gpfa | egrep -vw 'exon' | gawk 'NR>5 {OFS=\"\\t\";if (\$3==\"Gene\") print \"\#\$\"; \$2=\"geneid_$species\"; if (substr(\$1,1,1)!=\"\#\")
 print }' |";
    while (<LOCID>) {
	$geneidall .= $_;
    }
    close LOCID;

    my $tempgeneidgffpreds = $results.$species.".geneid.predictions.gff";
    open FOUT, ">$tempgeneidgffpreds";
    print FOUT "$geneidall";
    close FOUT;
    
    unlink $gff2psplots;
    
    #print STDERR "Plotting of predictions using gff2ps\n";

        open(LOCI,"<$gplen");
   	print STDERR "\nCreating gff2ps plots for $species\n\n";
  	 	while (<LOCI>) {
		    my ($gene_id,$genelength)=split;
		    	print STDERR "\n$gplen $gpgff $tempjkf_geneid $gene_id $genelength\n";
		    `egrep -w '$gene_id' $gpgff > $TMP/$gene_id.gff`;
		    `egrep -w '$gene_id' $tempgeneidgffpreds >> $TMP/$gene_id.gff`;
		    if ($jacknifevalidate){
		    `egrep -w '$gene_id' $tempjkf_geneid >> $TMP/$gene_id.gff`;
		}
		    if (!$contigopt){
		       `gff2ps -v -p -- $TMP/$gene_id.gff > $plotsdir/$species.${gene_id}.ps`;
		    print STDERR "#";}
		    elsif ($contigopt) {
			my $nucleotidesperline = 10000;
			`gff2ps -v -p -N $nucleotidesperline -C $path/.gff2psrcNEW -- $TMP/$gene_id.gff > $plotsdir/$species.ps`;
			`ps2pdf $plotsdir/$species.ps $plotsdir/$species.pdf`;
			`rm $plotsdir/$species.ps`;
			print STDERR "#";
		    }
    	}		       
    	close LOCI;


}


sub TblToFasta{

my($tbl,$faout) =  @_;

open(IN,"<$tbl");
open(FOUT,">$faout");
while (<IN>) {
	chomp $_;
	my ($n,$s) = split(/\s+/, $_);
	my ($i,$e) = (1, length($s));
	print FOUT ">$n\n";
	while ($i<=$e) {
		print FOUT substr($s,$i-1,60)."\n"; $i += 60;
	};
};

close IN;
close FOUT;

return $faout;

}

sub TblToFastaFile{

my($dir,$tbl) =  @_;

open(IN,"<$tbl");
while (<IN>) {
	chomp $_;
	my ($n,$s) = split(/\s+/o, $_); 
	open(FOUT,">${dir}"."$n");
	#print STDERR "${dir}"."$n";
	my ($i,$e) = (1, length($s));
	print FOUT ">$n\n";
	print STDERR "#";
	while ($i<=$e) {
		print FOUT substr($s,$i-1,60)."\n"; $i += 60;
	};
};


close IN;
close FOUT;

}


sub FastaToTbl{

my ($fa,$tblout) =  @_;

open(IN,"<$fa");
open(TOUT,">$tblout");
#print STDERR "$fa INSIDE LOOP\n\n";
#print STDERR "$tblout INSIDE LOOP\n\n";
 my $count = 0;
  while(<IN>){
    chomp;
    $_=~ s/\|//;
    if ( $_=~/\>(\S+)/ ){
      print TOUT "\n" if $count > 0;
      print TOUT $1."\t";
      $count++;
      }
    else{
    print TOUT $_;
      }
  }
print TOUT "\n";

close IN;
close TOUT;

return $tblout;

}


sub Translate{

my ($geneticcode,$cds,$outprot) = @_;
print STDERR "$geneticcode in loop\n";
my $frame = 0;

if (!open(FILEIN,"<$geneticcode")) {
  print "Translate: impossible to open genetic.code\n";
  exit(1);
}

my %gencodeh = ();

while (<FILEIN>) {

my $line = $_;

my ($aa,$codon) = split (/\s+/,$line);
#print STDERR "$codon\n";

$gencodeh{$codon} = $aa;

#print STDERR "HERE: $gencodeh{$codon}\n";

}

close FILEIN;


if (!open(CDSIN,"<$cds")) {
  print "Translate: impossible to open $cds\n";
  exit(1);
}

open POUT, ">$outprot";

while (<CDSIN>) {

my $line = $_;

my ($name,$seq) = split(/\s+/,$line);

my $lseq = length ($seq);

print POUT "$name ";

for (my $i=$frame;$i<$lseq-1;$i+=3) {
	  my $triplet = substr($seq,$i,3);
          if ($gencodeh{$triplet}){
	      print POUT "$gencodeh{$triplet}";
	  }else
	      {print POUT "X"; }

	} #for
	 print POUT "\n";
}

close CDSIN;
close POUT;

return $outprot;

}


#CONVERT GFF2 TO GENEID GFF
sub generalGFFtoGFFgeneid {

my($gff,$species,$type,$output) =  @_;
my %G;
my @G=();

my $geneidgff = ${output}.$species.${type}.".geneid_gff";

open(GFF,"<$gff");
open(GFFOUT,">$geneidgff");
while (<GFF>) {
  my ($c, @f, $id);
  $c = ":";
  $_=~ s/\|//;
  chomp;
  @f = split /\s+/o, $_;
  $id = $f[8]; #seq name i.e. 7000000188934730
  (exists($G{$id})) || do {
      $c = "#";
      $G{$id} = [ @f[8,0,6], 0, @f[3,4], [] ]; # [7000000188934730 1.4 - 0 46549	46680 ]
  };
  push @{ $G{$id}[6] }, [ @f[3,4] ];
  $G{$id}[3]++;
  $G{$id}[4] > $f[3] && ($G{$id}[4] = $f[3]);
  $G{$id}[5] < $f[4] && ($G{$id}[5] = $f[4]);
 };

@G = sort {  $a->[1] cmp $b->[1]
          || $a->[4] <=> $b->[4]
          || $a->[5] <=> $b->[5] }
      map { $G{$_} } keys %G;

 foreach my $GN (@G) {
    my ($id, $ctg, $str, $nex, $go, $ge,
        $coords, @coords, $ce, $CE, $cur, $feat, $c);
    ($id, $ctg, $str, $nex, $go, $ge, $coords) = @$GN;
    # print STDERR Data::Dumper->Dump([ \@$coords ], [ qw/ *coords / ])."\n";
    @coords =  map { $_->[0], $_->[1] }
              sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1]}
               map { $_ } @$coords;
    # print STDERR Data::Dumper->Dump([ \@coords ], [ qw/ *Coords / ])."\n";
    #print "# $id $ctg $str $nex $go $ge\n";
    $ce = 0; $CE = $nex * 2; $c = 1;
    while ($ce < $CE) {
        # $cur = ($str eq '-' ? $CE - $ce - 2 : $ce);
        if ($nex == 1) {
            $feat = "Single";
        } elsif ($c == 1) {
            $feat = $str eq '-' ? "Terminal" : "First";
        } elsif ($c == $nex) {
            $feat = $str eq '-' ? "First" : "Terminal";
        } else {
            $feat = "Internal";
        };
        print GFFOUT join("\t",$ctg,"$species",$feat,$coords[$ce],$coords[($ce + 1)],
                        ".",$str,".",$id)."\n";
        $ce += 2; $c++;
    };
};


my $geneidgffsorted = "";
open LOCID2, "sort -s -k8,9 -k4,5n $geneidgff |";
  	      while (<LOCID2>) {
  		  $geneidgffsorted .= $_;
  	      }
  	    close LOCID2;

	  
 	  ####
 	    my $tempgeneidgffsorted =  ${output}.$species.$type.".geneid.gff_sorted";   
 	    open FOUT2, ">$tempgeneidgffsorted";
 	    print FOUT2 "$geneidgffsorted";
 	    close FOUT2;


return $tempgeneidgffsorted;

close GFF;
close GFFOUT;
#exit(0);




}

sub go_to_die() { 
   (warn "@_\n") && &clean_tmp();
   &clean_ext(); #unless $ext_flg;
   exit(1); 
}


###EVAL OPTIMIZATION SORTING FUNCTION

sub clean_ext() {
    # Obtaining the list of extended files
    my @files = (qw( "" ));

   # Unlinking the temporary files if they exist
	foreach my $file (@files){
		unlink $file if (-e $file);
	}    
   # rmdir "geneid_params" if (-e "geneid_params");

}

sub clean_tmp() {
    # Obtaining the list of temporary files
    opendir(DIR,"$TMP");
    my @files = map { "$TMP/$_" } grep { /^$TMPROOT/ } readdir(DIR);
    closedir(DIR);

	foreach my $file (@files){
		unlink $file if (-e $file);
	}    
}


sub sorteval {

    $b->[7] <=> $a->[7]
      ||
	$b->[10] <=> $a->[10]
	  ||
	    $b->[4] <=> $a->[4]
	      ||
		$a->[11] <=> $b->[11]
		  ||
		    $a->[12] <=> $b->[12]

}


sub sortevalbranch {

    $b->[9] <=> $a->[9]
      ||
	$b->[12] <=> $a->[12]
	  ||
	    $b->[6] <=> $a->[6]
	      ||
		$a->[13] <=> $b->[13]
		  ||
		    $a->[14] <=> $b->[14]

}
