#!/usr/bin/gawk -f
#USAGE= $BIN/preparetrimatrixstart4parameter.awk 4 5 6 7 $DATA/set2/sites/starts.geneid.dimatrix
BEGIN {
    st=ARGV[1];
    nd=ARGV[2];
    rd=ARGV[3];
    th=ARGV[4];   
    ARGV[1]=ARGV[2]=ARGV[3]=ARGV[4]="";
}
{
  if (($1==st && $2!~/A$/) || ($1==nd && $2!="AT") || ($1==rd  && $2!="TG")  || ($1==th && $2!~/^G/)){
         print $1, $2, -9999;
    } else if (($1==st && $2~/A$/) || ($1==nd && $2=="AT") || ($1==rd  && $2=="TG"))
    print $1, $2, 0;
  else
print;
}
