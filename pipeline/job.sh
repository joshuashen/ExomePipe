#!/bin/sh
#$ -cwd
#$ -l mem=15G,time=16::

rm *.output

./run_pipeline.sh -I /ifs/scratch/c2b2/ip_lab/yshen/Yale/$1 -R /ifs/data/c2b2/ip_lab/shares/DATA/Sequencing/resources/bcm_hg18.fasta -E /ifs/scratch/c2b2/ip_lab/aps2157/ExomePipe/pipeline/data/old_samples/captureexoncoordinatesonly.txt -D /ifs/scratch/c2b2/ip_lab/aps2157/ExomePipe/pipeline/data/old_samples/dbsnp_129_hg18.ordered.rod > pipeline.output

