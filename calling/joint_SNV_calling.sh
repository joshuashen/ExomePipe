#!/bin/bash
#$ -cwd

HEAP=4
### TEMP="/ifs/scratch/c2b2/af_lab/saec/temp/"
INP=""
heap=20
targeted=1


USAGE="Usage: $0 -i <list of bam files> -m <heap> -s <global setting> [-t 1/0]"
ERRORMESSAGE="#### ERROR"
ERRORMESSAGE1="The following error has occurred"


while getopts i:m:o:t:s:h opt
  do 
  case "$opt" in
      i) INP="$OPTARG";;
      m) MEM="$OPTARG";;
      s) setting="$OPTARG";;  # global config
      t) targeted="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $INP == ""  || $setting == "" ]]
then
    echo $USAGE
    exit 1
fi

if [[ ! $MEM == "" ]]; then
    heap=$MEM
fi


. $setting

echo $ExonFile

outprefix=$INP."SNV_joint"

if [[ $REFTYPE != 'hg' ]]; then
    chrprefix=''
fi

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
  mkdir -p "temp/chr"$i
  
  out=${outprefix}"_chr"$i".sh"
  chrtarget=$out".list"

  echo $out
  echo $chrtarget

  egrep "^${i}" -w $ExonFile | awk '{print $1":"$2"-"$3}' > $chrtarget

  echo '#!/bin/bash'  > $out
  echo '#$ -cwd' >> $out
  
  
  
  cmd="java -Xmx${heap}g -Djava.io.tmpdir=temp/chr${i}/  -jar $GATKJAR -T UnifiedGenotyper  -R $REF  -D $DBSNP  -nt 2 -o ${outprefix}_chr${i}.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov 50 -L $chrtarget"
#   echo $cmd >> $out
  for bam in `cat $INP`
    do
    cmd=$cmd" -I "${bam}
  done
  echo $cmd >> $out

done
