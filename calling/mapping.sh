#!/bin/bash
#$ -cwd

# default values
ref="/ifs/data/c2b2/ip_lab/shares/DATA/Sequencing/resources/human_g1k_v37.fasta"
maxgaps=1
qualtrim=5
platform="illumina"
threads=2
sampleName=""
bwa=`which bwa`
samtools=`which samtools`
output=""
bwaversion=`$bwa 2>&1 | grep Version | awk '{print $2""$3}'`
sortmem=1000000000  # mem allocated for samtools sort  

### Note: bwa sample uses about 3.5G RAM,  if this job is submitted with 8G max RAM, the max $sortmem should be about 3G. 

USAGE="Usage: $0 -r ref -i foo_1.fastq [ -p foo_2.fastq ] [ -g maxgaps] [ -q qualtrim ] [ -n sampleName] [ -f platform] [ -s global_setting ]"

while getopts r:i:p:g:q:n:s:f:m:o:h opt
  do      
  case "$opt" in
      r) ref="$OPTARG";;
      i) fastq1="$OPTARG";;
      p) fastq2="$OPTARG";;
      m) sortmem="$OPTARG";;
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
    
# view -bS -o #{output} -

    $bwa sampe -p $platform -i $readgroup -l $readgroup -m $sampleName $ref $fastq1.sai $fastq2.sai $fastq1 $fastq2 | $samtools view -bS - | $samtools sort -m $sortmem  -  $output
    # echo $cmd
    rm -f $fastq1.sai $fastq2.sai
   
else
    
    $bwa samse -p $platform -i $readgroup -l $readgroup -m $sampleName $ref $fastq1.sai $fastq1 | $samtools view -bS - | $samtools sort -m $sortmem  -  $output
    rm -f $fastq1.sai
fi

# fix a bug in samtools sort
$samtools view -H $output.bam | sed 's/SO\:unsorted/SO:coordinate/' > $output.bam.header

# add mapping tool to header

echo -e "@PG\tID:$readgroup\tVN:$bwaversion\tCL:$bwa" >> $output.bam.header
# alternative for older version of echo:
# echo  "@PG"$'\t'"ID:$readgroup"$'\t'"VN:$bwaversion"$'\t'"CL:$bwa" >> $output.bam.header

$samtools reheader $output.bam.header $output.bam > $output.bam.temp 
mv $output.bam.temp $output.bam

# index
$samtools index $output.bam
