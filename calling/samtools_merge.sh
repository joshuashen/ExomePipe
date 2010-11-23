#!/bin/bash
#$ -cwd


GLOBAL=$1
prefix=$2

.  $GLOBAL

${SAMTOOLS} merge ${prefix}.realigned.bam ${prefix}.*.fixed.bam
${SAMTOOLS} index $prefix.realigned.bam 

