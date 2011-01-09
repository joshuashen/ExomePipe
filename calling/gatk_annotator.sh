#!/bin/bash
#$ -cwd

heap=4

while getopts v:g:t:m:h opt
  do  
  case "$opt" in
      v) vcf="$OPTARG";;
      g) GLOBAL="$OPTARG";;
      t) TEMP="$OPTARG";;
      m) MEM="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $vcf == "" || $GLOBAL == "" || $TEMP == "" ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL

if [[ ! $MEM == "" ]]; then
    heap=$MEM
fi

echo "heap size: ${heap}g"

targets=$vcf".targets.list"

awk '{print $1":"$2"-"$3}' $ExonFile > $targets

JAVA="java -Xmx${heap}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}

echo $GATK

$GATK \
    -T GenomicAnnotator \
    -R $REF \
    -B:variant,vcf $vcf \
    -L $targets \
    -B:refseq,AnnotatorInputTable $AnnotationTable \
    -m \
    -o $vcf.annotated \
    -BTI variant \

