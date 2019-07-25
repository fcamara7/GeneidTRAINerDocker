#!/usr/bin/gawk -f

# converts geneid gff output into DAS mysql compliance data

BEGIN{
	OFS=" ";
    }

{

 
    if (NR==1) {tot=split($0,array," ");seq=$1;add=$2;for(coord=4;coord<=tot;coord++){s=s""array[coord]" "};tot=0;print species,leng,$3,s;s="";oldadd=0;}
    if ($1!=seq) {tot=split($0,array," ");seq=$1;for(coord=4;coord<=tot;coord++){val=array[coord]+add;s=s""val" "};print species,leng,$3,s;s="";add=(add+$2);tot=0;}
   


 
}

