#!/bin/sh
#$ -S /bin/sh
#$ -cwd

# Filter SNPs in indel regions
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


### python $STING/python/makeIndelMask.py $INP.indels.raw.bed 10 $INP.indels.mask.bed

##if [[ $? != 0 ]]
#then
#        echo "Variant filtration: MakeIndelMask FAILED"
#        exit 1
# fi

$GATK \
 -T VariantFiltration \
 -R $REF \
 -o $INP.filtered.vcf \
 -B:variant,VCF $INP \
 -B:mask,VCF ${INDELVCF} \
 --maskName InDel \
 --clusterWindowSize 10 \
 --filterExpression "QUAL < 30.0 || QD < 5.0 || HRun > 5 || SB > -0.10" \
 --filterName "StandardFilters" \
 --filterExpression "MQ0 >= 4 && (MQ0 / (1.0 * DP)) > 0.1" \
 --filterName "HARD_TO_VALIDATE" 

### java -jar /path/to/dist/GenomeAnalysisTK.jar \
#  -T VariantFiltration \
#  -R /seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta \
#  -o /path/to/output.vcf \
#  -B:variant,VCF /path/to/input.vcf \
#  -B:mask,VCF /path/to/indels.calls.vcf \
#  --maskName InDel \
#  --clusterWindowSize 10 \
#  --filterExpression "AB < 0.2 || MQ0 > 50" \
#   --filterName "Nov09filters"