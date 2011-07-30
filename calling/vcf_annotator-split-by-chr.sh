#!/bin/bash
#$ -cwd

HEAP=2900

njobs=4  
## number of processes -- this module does not support multiple threads, so we have to 
## use background process to speed up things. 
## note: each gatk annotator instance uses about 3G RAM. 4 processes requires 12G. 

USAGE="Usage: $0 -v vcf -g global_settings -b bam.list [ -o output] [ -m heap_size] [-n njobs]"

while getopts v:b:o:g:n:m:h opt
  do  
  case "$opt" in
      v) input="$OPTARG";;
      b) bam="$OPTARG";;
      g) settings="$OPTARG";;
      o) output="$OPTARG";;
      n) njobs="$OPTARG";;
      m) HEAP="$OPTARG";;
      h)      echo $USAGE
	  exit 1;;
  esac
done

if [[ $input == "" || $settings == "" || $bam == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $settings 

temp=$input"_annotator_temp"
mkdir -p $temp

if [[ $output == "" ]] 
    then output=$input".moreinfo.vcf"
fi


JAVA="java -Xmx${HEAP}m -Djava.io.tmpdir="${temp}
GATK="$JAVA -jar "${GATKJAR}

uname -a 

w | grep load

cmdShared="$GATK  -T VariantAnnotator  -l INFO -R $REF -I $bam  -all  -B:variant,VCF $input  -B:dbsnp,VCF $DBSNPVCF  -B:compHapMap,VCF $HapMapV3VCF  -B:compdbSNP132,VCF $DBSNP132"

#    -B:comp1KG_CEU,VCF 1000GenomesCalls.CEU.vcf \

#    -A AlleleBalance \
#    -A DepthOfCoverage \
#    -A BaseQualityRankSumTest \
#    -A HomopolymerRun \
#    -A MappingQualityRankSumTest \
#    -A MappingQualityZero \
#    -A QualByDepth \
#    -A RMSMappingQuality \
#    -A SpanningDeletions \
#    -A HaplotypeScore  \


## split target list
target=$temp"/all-targets.list"

awk '{print $1":"$2"-"$3}' $ExonFile > $target
    
num=`wc -l $target | awk '{print \$1}'`

let nslice=$num/$njobs+1

total=0
for (( j=1; j<=$njobs; j++ ))  #  
  do

  chrtarget=$temp"/all-targets.slice."$j".list"

  let total=$total+$nslice
  
  if [[ $j -eq $njobs ]]
      then
      let njobs=$njobs-1
      let nslice=$num-$nslice*$njobs
      let njobs=$njobs+1
  fi
  head -${total} $target | tail -${nslice} > $chrtarget


  cmd=$cmdShared"  -L $chrtarget -o $output.$j.vcf"
  echo $cmd
  $cmd &  ## run it in background
done

wait

## combine

head -500 $output.1.vcf | egrep -w "^#" > $output

for (( j=1; j<=$njobs; j++ ))
  do
  egrep -w "^#" -v $output.$j.vcf >> $output
  rm -f $output.$j.vcf
done

rm -rf $temp
