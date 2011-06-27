#!/bin/sh
#$ -S /bin/sh
#$ -cwd

## filter based on allele balance -- which is good for most sites except highly polymorphic ones

HEAP=2
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
    -o $INP.filtere-AB.vcf \
    -B:variant,VCF $INP \
    --filterExpression "AB < 0.2" \
    --filterName "AlleleBalance"

rm -rf $TEMP

### java -jar /path/to/dist/GenomeAnalysisTK.jar \
#  -T VariantFiltration \
#  -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta \
#  -o /path/to/output.vcf \
#  -B:variant,VCF /path/to/input.vcf \
#  -B:mask,VCF /path/to/indels.calls.vcf \
#  --maskName InDel \
#  --clusterWindowSize 10 \
# 
#   --filterName "Nov09filters"