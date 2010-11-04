#!/bin/sh
#$ -cwd
#$ -l mem=15G,time=16::

#./run_pipeline.sh -I ./data/TCGA-13-1482-10A-01W-0549-09_IlluminaGA-DNASeq_exome.bam.sorted.bam -R ./data/wu_build36.fasta -L "1:20138-1268026" -D ./data/dbsnp_129_b36.ordered.rod -S /ifs/data/c2b2/ip_lab/shares/DATA/Sequencing/resources/refseq.autosome.rod -P illumina > pipeline.output

rm *.output

./run_pipeline.sh -I ./data/TCGA-13-1482-10A-01W-0549-09_IlluminaGA-DNASeq_exome.bam.sorted.bam -R ./data/wu_build36.fasta -E data/captureexoncoordinatesonly.txt -D ./data/dbsnp_129_b36.ordered.rod > pipeline.output

