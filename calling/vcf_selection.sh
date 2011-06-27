#!/bin/bash
#$ -cwd

### filter VCF based on:
# 1. calling rate (AN / 2 / N_samples)
# 2. ratio of zero-mapping-qual reads ( MQ0 / DP) or RMS of MQ 
# 3. allele frequency (AF or AC / 2 / N_samples)
#  etc

#default values
missing=0.1
mq=20
zmq=0.2 
freq=-1  # default no filter
ac=0  
maxfreq=1 # default no filter

HEAP=4

USAGE="Usage: $0 -v vcf -g global_settings [ -m missingness -q MQ -c min_allele_count -f min_Allele_freq -z zero_mq_ratio -F maxfreq ]"

while getopts v:m:q:z:f:F:s:h o
  do  
  case "$o" in
      v) vcf="$OPTARG";;
      m) missing="$OPTARG";;
      c) ac="$OPTARG";;
      q) mq="$OPTARG";;
      z) zmq="$OPTARG";;
      f) freq="$OPTARG";;
      F) maxfreq="$OPTARG";;
      s) samples="$OPTARG";;
      g) settings="$OPTARG";;
      h)      echo $USAGE
	  exit 1;;
  esac
done

if [[ $vcf == "" || $settings == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $settings 

nsample=`head -100 $vcf | grep "^#CHROM" | sed 's/\t/\n/g' | wc -l`
let nsample=$nsample-9

temp=$vcf"_filter_temp"
mkdir -p $temp
JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${temp}
GATK="$JAVA -jar "${GATKJAR}

# selection
selections=""
if [[ $samples != "" ]]; then
    selections=$selections" -sn $samples "
fi

let an=$nsample*2*(1-$missing)
selections=$selections" -select \"AN > $an \" -select \" AF > $freq \" -select \" AC > $ac \" -select \" AF < $maxfreq\" -select \"MQ > $mq  \"  "


java -jar GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R $REF \
    -B:variant,VCF $vcf \
    $selections \
    -env \
    -ef \
    -o $vcf.filtered

rm -rf $temp
