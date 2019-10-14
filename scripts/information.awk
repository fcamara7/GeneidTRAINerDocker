#! /usr/bin/gawk -f
# USAGE information.awk bakground_P P 

BEGIN {
    while (getline<ARGV[1]>0) 
	BG[$1,$2]=$4;
    ARGV[1]="";
    getline<ARGV[2];
    k=length($1);

}

{
    #print $1,$2,$4==1 ? 1 : $4*log($4/BG[$1,$2]); 
    if ($4==1) {print $1,$2,"1"}
    else if ($4==0) {print $1,$2,"0"}
    else if ($4!=0 && $4!=1) {print $1,$2,$4*log($4/BG[$1,$2])}
    #   print $1,$2,BG[$1,$2]==0 ? 0 : $4*log($4/BG[$1,$2]); 
    # print $1,$2,BG[$1,$2]==0 ? 0 : $4*log($4/BG[$1,$2])/log(2); 
   # c+=BG[$1,$2] == 0 ? 0 : $4*log($4/BG[$1,$2])/log(2);
   
    if ($4==1) {c+=1}
    else if ($4==0) {c+=0} 
    else if ($4!=0 && $4!=1) {c+=$4*log($4/BG[$1,$2])} 
    # c+=($4 == 0 || $4 == 1) ? 1 : $4*log($4/BG[$1,$2]);
    #  c+=($4 == 0 && $4 != 1) ? 0 : $4*log($4/BG[$1,$2]);
    #  c+=($4 == 1 && $4 != 0) ? 1 : $4*log($4/BG[$1,$2]);
    #c+=(BG[$1,$2] == 0 || $4 == 0) ? 0 : $4*log($4/BG[$1,$2]);
    #c+=(BG[$1,$2] == 0 || $4 == 0) ? 0 : $4*log($4/BG[$1,$2])/log(2);
    if (NR%(4^k) == 0) {
	#print $2,($4 == 0 || $4 == 1) ? c/2 : c;
	print $2,c;
 ;
	c=0;
    }
}
