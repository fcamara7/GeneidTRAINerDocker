#!/usr/bin/gawk -f
# Getkmatrix.awk
# gpf, imim, January 2000
# usage=gawk -f Getkmatrix.awk order number_of_nucleotides_per_sequence cds_seqfile (.tbl)

BEGIN {
  PCOUNT=0.25;

  k=ARGV[1];          #order of the markov chain
  num=ARGV[2];
  ARGV[1]=ARGV[2]="";

  # set pseudocounts 

  alpha[1]="A";
  alpha[2]="C";
  alpha[3]="G";
  alpha[4]="T";

  # for all k-tuples 
  for (i=1;i<=4;i++) 
     nx(1,k,alpha[i],N0,4*PCOUNT,num-k);


  # for all k+1-tuples
  for (i=1;i<=4;i++) 
     nx(1,k+1,alpha[i],N,PCOUNT,num-k);


#  sizek1=4^(k+1);

}
{
  sequence=toupper($2);
  if (sequence !~ /[^ACGT]/) { # consider only standard acgt

    lseq=length(sequence); 
    L+=(lseq-k);

    for (i=1;i<=lseq-k;i++) {
      ktuple=substr(sequence,i,k);
      ktuple1=substr(sequence,i,k+1);
      N0[i,ktuple]++; # ktuple frequence at position i
      N[i,ktuple1]++; # ktuple1 frequence at position i
      
    }
    total_seq++
  }
}
END {
  # get number of k-tuples observed in each frame
  L0=L+(3*sizek);

  for (t in N) {        # transition probabilities
    split (t,x, SUBSEP);
    pos=x[1];
    tk=substr(x[2],1,k);
    tk1=x[2];
    print x[1],x[2],N[t]/N0[pos,tk];
  }
}

function nx(l,len,s,Mat,p,num,   i,pos) { 
  if (l==len) {
    for (pos=1;pos<=num;pos++)
      Mat[pos,s]=p;
  }
  else {
     l++;
     for (i=1;i<=4;i++) 
       nx(l,len,s alpha[i], Mat,p,num);
   }
}
  
    








