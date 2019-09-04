# GeneidTRAINerDocker
Docker container containing the perl script that we use to automatically train the ab initio program geneid. It also creates a set of files used to evaluate geneid and other gene prediction programs that used in our #in house# protein-coding genome annotation pipeline (which will automated using NextFlow in the next few months)  

To start the docker container GeneidTRAINerDocker type the following:

docker build -t  geneidtrainerdocker .

In order to run it (the example below uses test data as input):

docker run -u $(id -u):$(id -g) -v $(pwd):/data -w /data geneidtrainerdocker -species M.cingulata -gff ./input/M.cingulata.cDNAs.450nt.complete.Uniprot98span.cds.4training4testing.gff2 -fastas ./input/M.cingulata.4training.fa -results ./output/ -sout stats.txt -branch no -reduced no 


