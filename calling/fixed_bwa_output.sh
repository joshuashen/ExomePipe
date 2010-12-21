#!/bin/bash
#$ -cwd

# do sort, rehead, and index

input=$1
header=$input.sort.bam.header

samtools sort -m 4000000000 $input $input.sort

samtools view -H $input.sort.bam | sed 's/SO\:unsorted/SO:coordinate/' > $header

samtools reheader $header $input.sort.bam > $input.sort.bam.temp

mv $input.sort.bam.temp $input.sort.bam

samtools index $input.sort.bam

### note: $input should be deleted later when you've confirmed sorting worked. 

