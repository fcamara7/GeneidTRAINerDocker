#USAGE= gawk -f aux.awk Background_P_file, Observed_P_file
BEGIN{
  while(getline<ARGV[1]>0) # read background probabilities
    BP[$1 $2]=$3;

  ARGV[1]="";
}
{
  if ($3!=0 && BP[$1 $2]!=0 ) print $1, $2, log($3/BP[$1 $2]);
  else print $1, $2, 0;
}
