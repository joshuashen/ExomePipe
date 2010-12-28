#!/bin/bash
#$ -cwd

# default values
ref="/ifs/data/c2b2/ip_lab/shares/DATA/Sequencing/resources/human_g1k_v37.fasta"
maxgaps=2
qualtrim=5
platform="illumina"
threads=2
sampleName=""
bwa=`which bwa`
samtools=`which samtools`
output=""
bwaversion=`$bwa 2>&1 | grep Version | awk '{print $2""$3}'`



USAGE="Usage: $0 -r ref -i foo_1.fastq [ -p foo_2.fastq ] [ -g maxgaps] [ -q qualtrim ] [ -n sampleName] [ -f platform] [ -s global_setting ]"

while getopts r:i:p:g:q:n:s:f:o:h opt
  do      
  case "$opt" in
      r) ref="$OPTARG";;
      i) fastq1="$OPTARG";;
      p) fastq2="$OPTARG";;
      g) maxgaps="$OPTARG";;
      q) qualtrim="$OPTARG";;
      n) sampleName="$OPTARG";;
      f) platform="$OPTARG";;
      o) output="$OPTARG";;
      s) setting="$OPTARG";;
      h)    echo $USAGE
          exit 1;;
  esac
done

if [[ $fastq1 == "" || ! -e $ref ]]; then
    echo $USAGE
    exit 1
fi

if [[ ! $setting == ""  ]]; then
    . $setting
fi




######### align step
cmd="$bwa aln -q $qualtrim -o $maxgaps -t  $threads  $ref  $fastq1 > $fastq1.sai"
echo $cmd
$bwa aln -q $qualtrim -o $maxgaps -t  $threads  $ref  $fastq1 > $fastq1.sai

readgroup=`basename $fastq1 | sed 's/.fastq$//'  | sed s'/.txt$//'`
if [[ $sampleName == "" ]]; then
    sampleName=$readgroup
fi

if [[ $output == "" ]]; then
    output=$fastq1.sorted
fi


if [[ ! $fastq2 == "" ]]; then  # paired-ends
    cmd="$bwa aln -q $qualtrim -o $maxgaps -t  $threads  $ref  $fastq2 > $fastq2.sai"
    echo $cmd
    $bwa aln -q $qualtrim -o $maxgaps -t  $threads  $ref  $fastq2 > $fastq2.sai
    
    $bwa sampe -p $platform -i $readgroup -l $readgroup -m $sampleName $ref $fastq1.sai $fastq2.sai $fastq1 $fastq2 | $samtools view -bS - | $samtools sort  -  $output
    # echo $cmd
    # $cmd
    
else
    
    $bwa sampe -p $platform -i $readgroup -l $readgroup -m $sampleName $ref $fastq1.sai $fastq1 | $samtools view -bS - | $samtools sort   -  $output
fi

# fix a bug in samtools sort
$samtools view -H $output.bam | sed 's/SO\:unsorted/SO:coordinate/' > $output.bam.header

# add mapping tool to header

echo "@PG\tID:$readgroup\tVN:$bwaversion\tCL:$bwa" >> $output.bam.header

$samtools reheader $output.bam.header $output.bam > $output.bam.temp 
mv $output.bam.temp $output.bam

# index
$samtools index $output.bam
