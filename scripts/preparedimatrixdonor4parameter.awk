#!/usr/bin/gawk -f
#USAGE= $BIN/preparedimatrixdonor4parameter.awk 3 4 5 $DATA/set2/sites/donor.geneid.dimatrix
BEGIN {
    st=ARGV[1];
    nd=ARGV[2];
    rd=ARGV[3];
    ARGV[1]=ARGV[2]=ARGV[3]="";
}
{
  if (($1==st && $2!~/G$/) || ($1==nd  && $2!="GT")  || ($1==rd && $2!~/^T/)){
         print $1, $2, -9999;
    } else if (($1==st && $2~/G$/) || ($1==nd  && $2=="GT"))
    print $1, $2, 0;
  else
print;
}
