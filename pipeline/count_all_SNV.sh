#!/bin/sh
#$ -cwd
# -l mem=5G,time=4::

GLOBAL="global_config.sh"

if [[ -e $GLOBAL ]]
then
        . $GLOBAL
else
        echo "Global config file not found. Exiting."
        exit 1
fi

USAGE="Usage: SNV_Count.sh -V <Variants file> -E <ExonFile> -I <Input Bam>"

while getopts V:E:I:h o
do      case "$o" in
        V)      VCF_FILE="$OPTARG";;
        E)      EXON_FILE="$OPTARG";;
        I)      INP="$OPTARG";;
        h)      echo $USAGE
                exit 1;;
        esac
done

if [[ $VCF_FILE == "" || $EXON_FILE == "" || $INP == "" ]]
then
        echo $USAGE
        exit 1
fi

DATAPATH=`dirname "$VCF_FILE"`
CONTIG_ORDER="`$SAMTOOLS idxstats $INP | grep -m 24 [0-9] | awk '{print $1}' | tr '\n' ' '`"

rm -f $DATAPATH/hap.list $DATAPATH/SNV.list
for i in $CONTIG_ORDER
do
	echo Test Printing
	./count_SNV.sh -V $VCF_FILE -E $EXON_FILE -L $i
	cat $DATAPATH/SNV.list.$i >> $DATAPATH/SNV.list
	rm -f $DATAPATH/SNV.list.$i $DATAPATH/hap.list.$i $DATAPATH/exon.list.$i
done


