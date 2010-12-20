#!/bin/sh
#$ -cwd

# samtools reheader <in.header.sam> <in.bam>

# deal with SO:unsorted after samtools sort

IN=$1

SAMTOOLS=`which samtools`

header=$IN.header


$SAMTOOLS view -H $IN | sed 's/SO\:unsorted/SO:coordinate/' > $header

$SAMTOOLS reheader $header $IN > $IN.temp

mv $IN.temp $IN

