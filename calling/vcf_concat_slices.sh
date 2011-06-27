#!/bin/bash
#$ -cwd

list=$1
out=$2

first=`head -1 $list`

head -200 $first | egrep "^#" > $out

##for (( i=1; i<=$n; i++ ))
for f in `cat $list`
  do 
  egrep -v "^#" $f >> $out
done
