#!/bin/bash
#$ -cwd

SAMTOOLS="/ifs/home/c2b2/ip_lab/yshen/usr/bin/samtools"

#### echo "sort bam"
####  ${SAMTOOLS} sort $INP $INP.sorted

${SAMTOOLS} index $1 
echo "$1 indexed"
