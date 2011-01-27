#!/bin/bash
#$ -cwd

### filter VCF based on:
# 1. calling rate (AN / 2 / N_samples)
# 2. ratio of zero-mapping-qual reads ( MQ0 / DP) or RMS of MQ 
# 3. allele frequency (AF or AC / 2 / N_samples)
#  etc

#default values
callrate=0.8
mq=20
zmq=0.4  # ((MQ0 / (1.0 * DP)) > $zmq 
freq=-1  # default no filter
ac=0  
maxfreq=1 # default no filter
# ab=0.95  # allele balance
qual=50.0 # min qual
clusterWinSize=10 #   --clusterWindowSize = 10
HRun=5 # homopolymer
qd=5.0 # qual over depth cutoff 
sb=-0.1 # strand bias
minDP=5 # average DP per sample, need to multiply by nsample

HEAP=2

USAGE="Usage: $0 -v vcf -g global_settings [ -w cluster_window_size -r callrate -q min_RMS_map-qual -c min_allele_count -f min_Allele_freq -z zero_mq_ratio -F maxfreq ] [ -i indelMask]"

while getopts v:r:q:z:f:F:s:w:Q:g:i:h o
  do  
  case "$o" in
      v) vcf="$OPTARG";;
      i) indelMask="$OPTARG";;
      r) callrate="$OPTARG";;
      c) ac="$OPTARG";;
      q) mq="$OPTARG";;
      z) zmq="$OPTARG";;
      f) freq="$OPTARG";;
      F) maxfreq="$OPTARG";;
      s) samples="$OPTARG";;
      w) clusterWinSize="$OPTARG";;
      Q) qd="$OPTARG";;
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
let minDP=$minDP*$nsample

temp=$vcf"_filter_temp"
mkdir -p $temp
JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${temp}
GATK="$JAVA -jar "${GATKJAR}

# selection
filter=""
if [[ $indelMask != "" ]]; then
    # indelMask is firstly a vcf, need to convert
    #echo $indelMask
    filter=$filter" -B:mask,VCF $indelMask  --maskName InDel "
fi

filter=$filter" --clusterWindowSize $clusterWinSize  "

uname -a 
echo  $GATK

let an=$nsample*2

## float point calculation in shell! 
an=`echo "scale=4; $callrate*$an" | bc`
#let an=$callrate*$an
echo $an

## AB filter is depreciated, set the threshold to be large (0.95); it's replaced by QD filter
echo $filter

# exit 

$GATK \
    -T VariantFiltration \
    -R $REF \
    -B:variant,VCF $vcf \
    $filter \
    -o $vcf.masked 

rm -rf $temp
