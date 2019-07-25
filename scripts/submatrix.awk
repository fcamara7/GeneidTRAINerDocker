#!/usr/bin/gawk -f

BEGIN {
    begin=ARGV[1];
    end=ARGV[2];
    ARGV[1]=ARGV[2]="";
}
($1 >= begin && $1 <= end ){
    $1=$1-begin+1;print $1,$2,$3;
}
