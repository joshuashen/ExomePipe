#!/bin/sh
#$ -cwd
#$ -l mem=5G,time=1::
# Findmem


# Load all samples from chromosome in a single $INP bam file

HEAP=4

INP=""
CHR=""
TEMP=""
GLOBAL=""
MEM=""

USAGE="Usage: $0 -I <Input bam file> -R <Reference Fasta>"

while getopts I:L:g:m:h o
do      case "$o" in
        I)      INP="$OPTARG";;
        L)      CHR="$OPTARG";;
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


TEMP=$INP"_temp-4-indel"

if [ ! -d $TEMP ]; 
    then    
    mkdir -p $TEMP
fi


JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}


cmd="$GATK -T IndelGenotyperV2 -R $REF -I $INP  -bed $INP.indels.brief.bed  -verbose $INP.indels.verbose.bed     -o $INP.indels.vcf  --refseq $REFSEQ  --maxNumberOfReads 50000 --window_size 300"
    
if [[ $CHR != "" ]]; then
    cmd=$cmd" -L $CHR"
fi

$cmd

if [[ $? == 0 ]]
    then
    echo "indel calling complete"
fi

rm -rf $TEMP

#-minCnt 2 \
#-minFraction 0.03 \
#-minConsensusFraction 0.6 \
#-mnr 1000000 
#-O $INP.indels.raw.bed \
#-o $INP.detalied.output.bed \
#-verbose \
