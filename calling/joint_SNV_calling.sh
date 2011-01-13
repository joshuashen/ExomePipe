#!/bin/bash
#$ -cwd

### TEMP="/ifs/scratch/c2b2/af_lab/saec/temp/"
heap=20  # default max heap size
targeted=1   # default: only call variants in targeted regions
nt=2 # default number of threads
dcov=300 # down sampling to dcov if depth is larger than dcov
mbq=20 # min base qual 
mmq=20 # min mapping qual

USAGE="Usage: $0 -i <list of bam files> -m <heap> -s <global setting> [-t 1/0] [ -n number_of_threads] [-d down_sampling] [ -b min_base_qual] [ -q min_mapping_qual]"
ERRORMESSAGE="#### ERROR"
ERRORMESSAGE1="The following error has occurred"


while getopts i:m:o:t:s:d:n:b:q:h opt
  do 
  case "$opt" in
      i) bamlist="$OPTARG";;
      m) MEM="$OPTARG";;
      s) setting="$OPTARG";;  # global config
      t) targeted="$OPTARG";;
      n) nthreads="$OPTARG";;
      d) dcov="$OPTARG";;
      b) mbq="$OPTARG";;
      q) mmq="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $bamlist == ""  || $setting == "" ]]
    then
    echo $USAGE
    exit 1
fi

if [[ ! $MEM == "" ]]; then
    heap=$MEM
fi

if [[ ! $nthreads == "" ]]; then
    nt=$nthreads
fi

. $setting

### get the name convention for chromosomes
chrString=`head -2 $REF | egrep "^>chr"`
if [[ $chrString != "" ]]; then
    chrprefix='chr'
else
    chrprefix=''
fi

temp=$bamlist"_SNVcall_dir"

# echo $ExonFile

outprefix=$bamlist."SNV_joint"


#target=$outprefix"_targets.all"

# awk '{print $1":"$2"-"$3}' $ExonFile > $target
    
for (( j=1; j<=24; j++ )) 
  do
  i=$j
  if [[ $i == "23" ]]; then
      i="X"
  fi
  if [[ $i == "24" ]]; then
      i="Y"
  fi
  echo "generating script for chr ${i}"
  mkdir -p "${temp}/chr"$i
  
  out=${outprefix}"_chr"$i".sh"
  chromosome=${chrprefix}${i}
  
  if [[ ! $targeted == "0" ]]; then  # only call variants within target intervals
      chrtarget=$temp"/targets."$chromosome".list"
      egrep "^${chromosome}" -w $ExonFile | awk '{print $1":"$2"-"$3}' > $chrtarget
  else  # call entire genome
      chrtarget=$chromosome
  fi
      
  
  echo '#!/bin/bash'  > $out
  echo '#$ -cwd' >> $out
  
  
  
  cmd="java -Xmx${heap}g -Djava.io.tmpdir=${temp}/chr${i}/  -jar $GATKJAR -T UnifiedGenotyper  -R $REF  -D $DBSNP  -nt ${nt} -o ${outprefix}_chr${i}.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov ${dcov} -mbq ${mbq}  -mmq ${mmq} -L $chrtarget -I $bamlist"


  
  echo $cmd >> $out
  
  let qmem=$heap+3
  qsub -l mem=${qmem}G,time=60:: -o $temp/log.$i.o -e $temp/log.$i.e -N joint.call.$bamlist.$i $out 
  # echo $qmem
done
