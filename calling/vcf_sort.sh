#!/bin/bash
#$ -cwd


vcf=$1
chrlist=$2

head -100 $vcf | egrep "^#" > $vcf.sorted

# testl=`head -100 $vcf | tail -1`

for f in `cat $chrlist` 
  do
  egrep -w "^$f" $vcf | sort -k2n >> $vcf.sorted
done

