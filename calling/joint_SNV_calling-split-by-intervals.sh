#!/bin/bash
#$ -cwd

nt=2 # default number of threads
dcov=300 # down sampling to dcov if depth is larger than dcov
mbq=20 # min base qual 
mmq=20 # min mapping qual
njobs=900

USAGE="Usage: $0 -i <list of bam files> -m <heap> -s <global setting> [ -n number_of_threads] [ -j number_of_qjobs] [-d down_sampling] [ -b min_base_qual] [ -q min_mapping_qual] [ -v total_mem ]"

while getopts i:m:o:s:d:n:j:b:q:v:h opt
  do 
  case "$opt" in
      i) bamlist="$OPTARG";;
      m) MEM="$OPTARG";;
      s) setting="$OPTARG";;  # global config
      n) nthreads="$OPTARG";;
      d) dcov="$OPTARG";;
      j) njobs="$OPTARG";;
      b) mbq="$OPTARG";;
      q) mmq="$OPTARG";;
      v) qmem="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $bamlist == ""  || $setting == "" ]]
    then
    echo $USAGE
    exit 1
fi

if [[ ! $nthreads == "" ]]; then
    nt=$nthreads
fi

let heap=$nt*4

if [[ ! $MEM == "" ]]; then
    heap=$MEM
fi

if [[ $qmem == "" ]]; then
    let qmem=$heap+3  # allocated RAM for each qjob in total
fi

. $setting

### get the name convention for chromosomes
# chrString=`head -2 $REF | egrep "^>chr"`
#if [[ $chrString != "" ]]; then
#    chrprefix='chr'
#else
#    chrprefix=''
# fi

bamlist=`readlink -e $bamlist`
temp=$bamlist"_SNVcall_dir"

mkdir -p $temp

# echo $ExonFile

target=$temp"/all-targets.list"

awk '{print $1":"$2"-"$3}' $ExonFile > $target
    
num=`wc -l $target | awk '{print \$1}'`

let nslice=$num/$njobs+1

bname=`basename $bamlist`
total=0
for (( j=1; j<=$njobs; j++ ))  #  
  do
  tempd=${temp}"/slices/"$j

  mkdir -p $tempd
  
  out=${temp}"/snv.slice."$j".sh"
#  chromosome=${chrprefix}${i}
  chrtarget=$out".targets.list" 

  let total=$total+$nslice

  if [[ $j -eq $njobs ]]
      then
      let njobs=$njobs-1
      let nslice=$num-$nslice*$njobs
      let njobs=$njobs+1
  fi
 
  head -${total} $target | tail -${nslice} > $chrtarget

  
  echo '#!/bin/bash'  > $out
  echo '#$ -cwd' >> $out
  echo 'uname -a' >> $out
  cmd="java -Xmx${heap}g -Djava.io.tmpdir=${tempd}  -jar $GATKJAR -T UnifiedGenotyper  -R $REF  -D $DBSNP  -nt ${nt} -o ${temp}/snv.slice.$j.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov ${dcov} -mbq ${mbq}  -mmq ${mmq} -L $chrtarget -I $bamlist"
  
  echo $cmd >> $out
  
  
  qsub -l mem=${qmem}G,time=60:: -o $temp/log.$j.o -e $temp/log.$j.e -N snp.$j.$bname $out 
  # echo $qmem
done
