#USAGE= gawk -f aux.awk Background_P_file, Observed_P_file
BEGIN{
  while(getline<ARGV[1]>0) # read backgrounb probabilities
    BP[$1 $2]=$4;

  ARGV[1]="";
}
{
  if ($4!=0 && BP[$1 $2]!=0 && $4!=1) print $1, $2, log($4/BP[$1 $2]);
  if ($4!=0 && BP[$1 $2]!=0 && $4==1) print $1, $2, 0;
  if ($4==0 && BP[$1 $2]!=0) print $1, $2, -9999;
#else print $1, $2, 0;
}
