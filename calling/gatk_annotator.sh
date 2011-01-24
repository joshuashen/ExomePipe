#!/bin/bash
#$ -cwd

heap=4

while getopts v:g:m:h opt
  do  
  case "$opt" in
      v) vcf="$OPTARG";;
      g) GLOBAL="$OPTARG";;
      m) MEM="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $vcf == "" || $GLOBAL == "" ]]
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

TEMP=$vcf"_anno-temp"
mkdir -p $TEMP

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

rm -rf $TEMP
