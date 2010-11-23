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
USAGE="Usage: $0 -I <Input bam file> -R <Reference fasta> -D <DBSNP file> -g <Global config>  [[-L #:#-#] | [-E <ExonFile>]]\n"


ERRORMESSAGE="#### ERROR"
ERRORMESSAGE1="The following error has occurred"

while getopts I:L:R:D:E:g:h o
do      case "$o" in
        I)      INP="$OPTARG";;
        L)      CHR="$OPTARG";;
        R)      REF="$OPTARG";;
        D)      DBSNP="$OPTARG";;
        E)      ExonFile="$OPTARG";;
        g)      GLOBAL="$OPTARG";;  # global config 
	h)	echo $USAGE
		exit 1;;
        esac
done

if [[ $INP == "" || $REF == "" || $DBSNP == "" || $GLOBAL == "" ]]
then
	echo $USAGE
	exit 1
fi

. $GLOBAL

echo `date`

echo "make directory"
OUTDIR=$INP"_pipe"
if [ ! -d $OUTDIR ]; then
  mkdir $OUTDIR
fi
echo $OUTDIR 

echo "Calling gatk_realign.scr step by step..."

if [[ $CHR == "" ]]
then
qsub -l mem=6G,time=24:: -t 1-24 -sync y -o $OUTDIR/realignment.stdout -e $OUTDIR/realignment.stderr ${BPATH}/gatk_realign.scr -g $GLOBAL -I $INP -R $REF -D $DBSNP -m 0
else
qsub -l mem=6G,time=24::  -sync y -o $OUTDIR/realignment.stdout -e $OUTDIR/realignment.stderr ${BPATH}/gatk_realign.scr -g $GLOBAL -I $INP -R $REF -D $DBSNP -m 0 -L $CHR
fi


## !!! need to check if realign X/Y chromosomes were complete (probably looking up the log file)

if [[ $? != 0 || `grep "$ERRORMESSAGE" $OUTDIR/realignment.std*` != "" ]]
    then
    echo "Realignment FAILED"
    exit 1
fi

rm -rf $OUTDIR/realignment.*.std*
cat $OUTDIR/*.forRealigner.intervals | bzip2 - > $OUTDIR/forRealigner.intervals.bz2

##  merge bam files
echo "merge bam files"

${SAMTOOLS} merge $INP.realigned.bam $OUTDIR/*.fixed.bam 
${SAMTOOLS} index $INP.realigned.bam

echo "clean house"
rm -rf ${INP}.*.fixed.bam 

# change ${INP}
$INP=${INP}.realigned.bam



echo "Execution of ./gatk_realign.scr complete"
echo

# rm -rf {$INP} 

# $INP=${INP}.sorted.bam



date
echo "Calling ./gatk_recalibrate.scr ..."
qsub -l mem=7G,time=16:: -sync y -o $OUTDIR/recalibrate.stdout -e $OUTDIR/recalibrate.stderr ${BPATH}/gatk_recalibrate.scr  -g $GLOBAL  -I $INP -R $REF

if [[ $? != 0 || `grep "$ERRORMESSAGE" $OUTDIR/recalibrate.std*` != "" ]]
then
	echo " Recalibration FAILED"
	exit 1
fi
echo "Execution of ./gatk_recalibrate.scr complete"
echo

$INP=${INP}.recalibrated.bam

date
echo "Calling ./gatk_depthofcoverage.scr ..."

qsub -l mem=5G,time=12:: -o $OUTDIR/depthofcoverage.stdout -e $OUTDIR/depthofcoverage.stderr ${BPATH}/gatk_depthofcoverage.scr  -g $GLOBAL -I $INP -R $REF -E "$ExonFile"

date
echo "Calling ./gatk_indelcall.scr ..."
qsub -l mem=5G,time=16::  -o $OUTDIR/indelcalling.stdout -e $OUTDIR/indelcalling.stderr ${BPATH}/gatk_indelcall.scr  -g $GLOBAL -I $INP -R $REF


if [[ $? != 0 || `grep "$ERRORMESSAGE" $OUTDIR/indelcalling.std*` != "" ]]
    then
    echo "Indel calling FAILED"
    exit 1
fi
echo "Execution of ./gatk_indelcall.scr complete"
echo



date
echo "Calling ./gatk_snpcall.scr ..."

qsub -l mem=5G,time=16:: -sync y -o $OUTDIR/snpcalling.stdout -e $OUTDIR/snpcalling.stderr ${BPATH}/gatk_snpcall.scr  -g $GLOBAL -I $INP -R $REF -D $DBSNP

if [[ $? != 0 || `grep "$ERRORMESSAGE" $OUTDIR/snpcalling.std*` != "" ]]
then
	echo "SNP calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_snpcall.scr complete"
echo



date



echo "Calling ./gatk_varianteval.scr ..."
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

