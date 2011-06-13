#!/bin/bash
#$ -cwd

# Filter indels
# Load all samples from chromosome in a single $INP vcf file


HEAP=4
INP=""
REF=""

USAGE="Usage: $0 -I <Input VCF file> -g "

while getopts I:t:m:g:h o
do      case "$o" in
        I)      INP="$OPTARG";;
        m)      MEM="$OPTARG";;
        g)      GLOBAL="$OPTARG";;
        h)      echo $USAGE
                exit 1;;
        esac
done

if [[ $INP == "" || $GLOBAL == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $GLOBAL


if [ ! $MEM == "" ]
    then
    HEAP=$MEM
fi

if [ $JOB_ID == "" ]; then
    JOB_ID="VarFilter"
fi

# JOB_ID is the qsub job ID
TEMP=$INP"_"$JOB_ID"/"

if [ ! -d $TEMP ]; 
    then    
    mkdir -p $TEMP
fi

JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}



$GATK \
    -T VariantFiltration \
    -R $REF \
    -o $INP.filtered.vcf \
    -B:variant,VCF $INP \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"  \
    --filterName "HARD_TO_VALIDATE" \
    --filterExpression "SB >= -1.0" \
    --filterName "StrandBiasFilter" \
    --filterExpression "QUAL < 10" \
    --filterName "QualFilter"