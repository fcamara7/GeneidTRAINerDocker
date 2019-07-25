#!/usr/bin/gawk -f
#USAGE= $BIN/preparedimatrixacceptor4parameter.awk 18 19 20 $DATA/set2/sites/acceptor.geneid.dimatrix
BEGIN {
    st=ARGV[1];
    nd=ARGV[2];
    rd=ARGV[3];
    ARGV[1]=ARGV[2]=ARGV[3]="";
}
{
  if (($1==st && $2!~/A$/) || ($1==nd  && $2!="AG")  || ($1==rd && $2!~/^G/)){
         print $1, $2, -9999;
    } else if (($1==st && $2~/A$/) || ($1==nd  && $2=="AG"))
    print $1, $2, 0;
  else
print;
}
