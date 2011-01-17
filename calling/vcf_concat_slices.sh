#!/bin/bash
#$ -cwd

n=$1
out=$2

head -200 snv.slice.1.raw.vcf | egrep "^#" > $out

for (( i=1; i<=$n; i++ ))
  do 
  egrep -v "^#" snv.slice.$i.raw.vcf >> $out
done
