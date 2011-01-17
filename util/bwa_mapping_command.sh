#!/bin/bash
#$ -cwd

## mapping paired reads; one run per batch.

date
# echo "align first batch"
/ifs/home/c2b2/ip_lab/yshen/usr/bin/bwa aln -t 2 /ifs/home/c2b2/ip_lab/yshen/data/ref_genomes/human_g1k_v37.fasta knome.1.fastq > knome.1.fastq.sai

date
echo "align second batch"
/ifs/home/c2b2/ip_lab/yshen/usr/bin/bwa aln -t 2 /ifs/home/c2b2/ip_lab/yshen/data/ref_genomes/human_g1k_v37.fasta knome.2.fastq > knome.2.fastq.sai

date
echo "pe"
/ifs/home/c2b2/ip_lab/yshen/usr/bin/bwa sampe  /ifs/home/c2b2/ip_lab/yshen/data/ref_genomes/human_g1k_v37.fasta knome.1.fastq.sai knome.2.fastq.sai knome.1.fastq knome.2.fastq   | bzip2 > knome.aligned.sam.bzip2

date
echo "done"


