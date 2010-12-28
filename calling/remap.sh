#!/bin/bash
#$ -cwd


# remap reads:
# 1. convert recalibrated BAMs to fastq 
# 2. bwa sampe 
# 3. realign 


# step 1. BAM to fastq

USAGE="Usage: $0 "

while getopts b:s:l:h opt
  do      
  case "$opt" in
      b) bam="$OPTARG";;
      s) setting="$OPTARG";;
      h)    echo $USAGE
          exit 1;;
  esac
done

if [[ $bam == "" || $setting == "" ]]; then
    echo $USAGE
    exit 1
fi

. $setting 

# get sampleName
sampleName=`$SAMTOOLS view -H $bam | grep '^@RG' |  sed 's/\s/\n/g' | grep '^SM:' | cut -f2 -d ':' | sort -u`


if [[ $sampleName == ""  ]]; then # no interested sample
    echo "no sample of interest"
    exit 1
else
    echo $sampleName
fi

# first, samtofastq
echo "${BPATH}/samtofastq.sh  -b $bam -p $sampleName -g $setting"
sh ${BPATH}/samtofastq.sh  -b $bam -p $sampleName -g $setting   ## setting is only useful to get programs, not $REF

# then, map 
echo "${BPATH}/mapping.sh -r $REF -i $sampleName"_1_sequence.fastq" -p $sampleName"_2_sequence.fastq" -n $sampleName -o $sampleName.sorted.bam"

sh ${BPATH}/mapping.sh -r $REF -i $sampleName"_1_sequence.fastq" -p $sampleName"_2_sequence.fastq" -n $sampleName -o $sampleName.sorted

bzip2 $sampleName"_1_sequence.fastq" $sampleName"_2_sequence.fastq"

OUTDIR=$sampleName".sorted.bam_pipe"
if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi  

for (( i=1; i<=24; i++ ))
  do 
  echo "qsub -N realign.$sampleName.$i -l mem=5G,time=48:: -o $OUTDIR/log.$i.realign.o -e $OUTDIR/log.$i.realign.e ${BPATH}/gatk_realign_atomic.scr -I $sampleName.sorted.bam  -g $newsetting -L $i "
  qsub -N realign.$sampleName.$i -l mem=5G,time=48:: -o $OUTDIR/log.$i.realign.o -e $OUTDIR/log.$i.realign.e ${BPATH}/gatk_realign_atomic.scr -I $sampleName.sorted.bam  -g $newsetting -L $i
done



