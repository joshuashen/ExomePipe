#!/bin/bash
#$ -cwd

## concatenate VCFs from different chromosomes of same sample(s)

# take the header of the first VCF 

vcflist=$1

head -1 $vcflist | xargs head -100 | grep ^#
 
for f in `cat $vcflist`
  do
  cat $f | grep -v ^#
done



# (zcat A.vcf.gz | head -100 | grep ^#; \ 
# zcat A.vcf.gz | grep -v ^#; \ 
# zcat B.vcf.gz | grep -v ^#; ) \ 
# | bgzip -c > out.vcf.gz