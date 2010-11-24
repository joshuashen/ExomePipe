#!/bin/sh
#$ -cwd

GLOBAL="global_config.sh"

if [[ -e $GLOBAL ]]
then
	. $GLOBAL
else
        echo "Global config file not found. Exiting."
        exit 1
fi

if [[ $1 == "" ]]
then
	echo "Please pass input bam file as parameter"
#	exit 1
fi

INPUT="/ifs/scratch/c2b2/ip_lab/yshen/Yale/RNAseq19082010/bamfiles/0/s_10242_1_sequence.txt.bam.sorted.bam"
DATAPATH=`dirname "$INPUT"`
rm -f $DATAPATH/*.output

#BPATH/run_pipeline.sh -I $1 -R $WU_REF -E $ExonList -D $WU_DBSNP -P $WU_Platform > $DATAPATH/pipeline.output 
$BPATH/run_pipeline.sh -I $INPUT -R $YALE_REF -E $ExonList -D $YALE_DBSNP -P $YALE_Platform > $DATAPATH/pipeline.output

