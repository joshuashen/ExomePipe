#!/bin/bash
#$ -cwd

HEAP=4

while getopts i:s:m:h o
  do 
  case "$o" in
      i) bam="$OPTARG";;
      s) setting="$OPTARG";;
      m) mem="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done


if [[ $bam == "" ||  $setting == "" ]]; then
    echo "Usage: $0 -i foo.bam -s setting [ -m heap ]"
    exit 0
fi

.  $setting

if [[ ! $MEM == "" ]]; then
    HEAP=$MEM
fi

temp=$bam".temp_dir/"

mkdir -p $temp

JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${temp}
MarkDup="$JAVA -jar ${PICARD}/MarkDuplicates.jar"

$MarkDup I=$bam O=$bam.temp METRICS_FILE=$bam.dup VALIDATION_STRINGENCY=SILENT

if [[ $? == 0 ]]
    then
    mv $bam.temp $bam
    rm $bam.bai
    $SAMTOOLS index $bam
    
    echo "dedup complete"
    rm -rf $temp
else
    echo "dedup failed"
    rm -rf $temp
    exit 1
fi
