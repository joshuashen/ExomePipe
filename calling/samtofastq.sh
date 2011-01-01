#!/bin/bash
#$ -cwd

heap=4 # 4G
prefix=""

USAGE="$0 -b foo.bam -p prefix_of_fastq -g global_config"

while getopts b:p:m:g:h o
  do      case "$o" in
        b)      bam="$OPTARG";;
        p)      pre="$OPTARG";;
        m)      mem="$OPTARG";;
        g)      global="$OPTARG";;  # global config
        h)      echo $USAGE
                exit 1;;
  esac
done

if [[ $bam == "" || $global == ""  ]]
    then
    echo $USAGE
    exit 1
fi

. $global

if [ ! $mem == "" ]
    then
    heap=$mem
fi

OUTDIR=$bam"_remap"
if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi


# JOB_ID is the qsub job ID
TEMP=$OUTDIR"/temp/sam2fastq/"

if [ ! -d $TEMP ]; then
    mkdir -p $TEMP
fi

if [ ! $pre == "" ]
    then
    prefix=$pre
else
    prefix=`basename $bam | sed 's/.bam$//'`
    prefix=$OUTDIR"/"$prefix
fi

# echo $prefix

JAVA="java -Xmx${heap}g -Djava.io.tmpdir="${TEMP}
SamToFastq="$JAVA -jar ${PICARD}/SamToFastq.jar" 


$SamToFastq INPUT=$bam FASTQ=$prefix"_1.fastq" SECOND_END_FASTQ=$prefix"_2.fastq"

echo "samtofastq complete"
rm -rf $TEMP


