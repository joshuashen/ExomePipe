#!/bin/bash
#$ -cwd

INP=""
CHR=""
REF=""
DBSNP=""
ExonFile=""
GLOBAL=""

## ExonFile is a tab-delimited file specifying the targeted regions ( such as the exome sequencing):
## format: chr start end
USAGE="Usage: $0 -I <Input bam file>  -g <Global config>\n"

while getopts I:g:h o
  do      
  case "$o" in
      I)    INP="$OPTARG";;
      g)    GLOBAL="$OPTARG";;  # global config 
      h)    echo $USAGE
	  exit 1;;
  esac
done

if [[ $INP == "" || $GLOBAL == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $GLOBAL


echo `date`

OUTDIR=$INP"_pipe"
if [ ! -d $OUTDIR ]; then
  mkdir $OUTDIR
fi

# JOB_ID is the qsub job ID
TEMP=$OUTDIR"/temp/"

if [ ! -d $TEMP ]; then
    mkdir -p $TEMP
fi

if [ ! -e $OUTDIR"/forRealigner.intervals.bz2" ]
    then
    cat $OUTDIR/*.forRealigner.intervals | bzip2 - > $OUTDIR/forRealigner.intervals.bz2
    rm -f $OUTDIR/*.forRealigner.intervals
fi

##  merge bam files

if [ ! -e $OUTDIR"/all.recalibrated.bam" ] # need recalib
    then
    if [ ! -e $OUTDIR"/all.realigned.bam"  ] # need merge
	then
	echo "merge bam files"
	${SAMTOOLS} merge $OUTDIR/all.realigned.bam $OUTDIR/*.fixed.bam  > $OUTDIR/log.merge 2>&1
	
	rm -f $OUTDIR/*.fixed.bam
	${SAMTOOLS} index $OUTDIR/all.realigned.bam  > $OUTDIR/log.index 2>&1
    fi
    date
    echo "recalibrate qual"
    
#  give it 8G heap 
    ${BPATH}/gatk_recalibrate.scr -m 8000  -g $GLOBAL -t $TEMP  -I $OUTDIR/all.realigned.bam  -o $OUTDIR/all.recalibrated.bam 
    if [[ $? == 0 ]]  # test if previous command completed
	then
	rm -f $OUTDIR/all.realigned.bam    
    else
	echo "recalibrate failed"
	exit 1
    fi
fi


date
cmd="qsub -N depth.$INP -l mem=9G,time=24:: -o $OUTDIR/log.depth.o -e $OUTDIR/log.depth.e ${BPATH}/gatk_depthofcoverage.scr  -g $GLOBAL -I $OUTDIR/all.recalibrated.bam -t $TEMP  -m 8000"
echo $cmd

$cmd


date

cmd="qsub -N indel.$INP -l mem=9G,time=24::  -o $OUTDIR/log.indel.o -e $OUTDIR/log.indel.e  ${BPATH}/gatk_indelcall.scr  -g $GLOBAL -I $OUTDIR/all.recalibrated.bam -t $TEMP -m 8000"
echo $cmd
$cmd

date
cmd="qsub -N snp.$INP -l mem=9G,time=24:: -o $OUTDIR/log.snp.o -e $OUTDIR/log.snp.e  ${BPATH}/gatk_snpcall.scr  -g $GLOBAL -I $OUTDIR/all.recalibrated.bam -t $TEMP -m 8000"
echo $cmd
$cmd


exit 0


date
echo "Evaluate SNPs"

### !!! Need indel calls !! 

# python $STING/python/makeIndelMask.py $INP.indels.r 10 $INP.indels.mask.bed

if [[ $? != 0 ]]
then
    echo "Variant Eval: MakeIndelMask FAILED"
    exit 1
fi

$GATK \
    -T VariantEval -R $REF \
    -B:eval,VCF $INP.snps.raw.vcf \
    -D $DBSNP \
    -E CountVariants \
    -noStandard \
    -o $INP.snps.raw.vcf.eval \ 
    -BTI eval


qsub -l mem=5G,time=1:: -sync y -o $OUTDIR/varianteval_prefiltering.stdout -e $OUTDIR/varianteval_prefiltering.stderr ${BPATH}/gatk_varianteval.scr  -g $GLOBAL -I $INP -R $REF -D $DBSNP
if [[ $? != 0 || `grep "$ERRORMESSAGE" $OUTDIR/varianteval_prefiltering.std*` != "" ]]
then
	echo "Prefilter variant eval calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_varianteval.scr complete"
echo



date
echo "Calling ./gatk_filter.scr ..."
qsub -l mem=5G,time=1:: -sync y -o $OUTDIR/variantfiltering.stdout -e $OUTDIR/variantfiltering.stderr ${BPATH}/gatk_filter.scr  -g $GLOBAL -I $INP -R $REF
if [[ $? != 0 || `grep "$ERRORMESSAGE" $OUTDIR/variantfiltering.std*` != "" ]]
then
	echo "Variant filtering FAILED"
	exit 1
fi
echo "Execution of ./gatk_filter.scr complete"
echo



date
echo "Calling ./gatk_varianteval.scr ..."
qsub -l mem=5G,time=1:: -sync y -o $OUTDIR/varianteval_postfiltering.stdout -e $OUTDIR/varianteval_postfiltering.stderr ${BPATH}/gatk_varianteval.scr  -g $GLOBAL -I $INP -R $REF -D $DBSNP
if [[ $? != 0 || `grep "$ERRORMESSAGE" $OUTDIR/varianteval_postfiltering.std*` != "" ]]
then
	echo "Postfilter variant eval calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_varianteval.scr complete"
echo
date

rm -fr $TEMP
