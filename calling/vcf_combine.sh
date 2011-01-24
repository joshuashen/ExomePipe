#!/bin/bash
#$ -cwd

USAGE="Usage: $0 -l List_of_VCFs -s settings [ -m heapSize]"
HEAP=4

while getopts l:s:m:h o
  do 
  case "$o" in
      l) list="$OPTARG";;
      m) MEM="$OPTARG";;
      s) GLOBAL="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $list == "" || $GLOBAL == "" ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL

if [[ $MEM != "" ]]
    then
    HEAP=$MEM
fi

inputline=""

# UNIQUIFY


for f in `cat $list`
do
  # sampleID=`echo $f | cut -f1 -d '_'`
  g=`echo $f | cut -f1 -d '.'`
  inputline=$inputline"-B:$g,VCF $f  "
done

# echo $inputline



temp=$list"_combine_temp"
mkdir -p $temp




JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${temp}
GATK="$JAVA -jar "${GATKJAR}

$GATK -T CombineVariants -R $REF  \
    -o $list"_combined.vcf"  \
    $inputline  \
    -variantMergeOptions UNION \
    -genotypeMergeOptions UNSORTED  \

rm -rf $temp