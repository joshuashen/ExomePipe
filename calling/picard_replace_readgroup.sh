#!/bin/bash
#$ -cwd

readgroup=""
sample=""
readlib=""
platform="illumina"

USAGE="Usage: sh _.sh -i input -o output -r readgroup -s sampleName -g global_setting\n"

while getopts i:g:o:r:s:l:h opt
  do     
  case "$opt" in
      i) input="$OPTARG";;
      o) output="$OPTARG";;
      r) readgroup="$OPTARG";;
      s) sample="$OPTARG";;
      l) readlib="$OPTARG";;
      g) GLOBAL="$OPTARG";;  # global config
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $readgroup == "" || $sample == ""  || $GLOBAL == "" ]]; then
    echo $USAGE
    exit 1
fi

if [[ $readlib == "" ]]; then
    readlib=$readgroup
fi

HEAP=1

. $GLOBAL


TEMP=$input"_replaceRG_temp"
if [ ! -d $TEMP ]; then
  mkdir $TEMP
fi


JAVA="java -Xmx${HEAP}m -Djava.io.tmpdir="${TEMP}

# FIXMATE="$JAVA -jar ${PICARD}/FixMateInformation.jar"

# $FIXMATE INPUT=$OUTDIR/$CHR.cleaned.bam OUTPUT=$OUTDIR/$CHR.fixed.bam SO=coordinate VALIDATION_STRINGENCY=SILENT

java -Xmx${HEAP}g   -Djava.io.tmpdir=${TEMP} \
    -jar ${PICARD}/AddOrReplaceReadGroups.jar  \
    RGLB=$readlib RGPL=$platform RGPU=$readlib RGSM=$sample RGID=$readgroup    \
    I=$input O=$output  \
    SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT

rm -rf $TEMP