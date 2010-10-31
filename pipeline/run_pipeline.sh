#!/bin/sh
#$ -cwd

INP=""
CHR=""
REF=""
DBSNP=""
ExonFile=""

USAGE="Usage: $0 -I <Input bam file> -R <Reference fasta> -D <DBSNP file> [[-L #:#-#] | [-E <ExonFile>]]"
ERRORMESSAGE="#### ERROR"
ERRORMESSAGE1="The following error has occurred"

while getopts I:L:R:D:E:h o
do      case "$o" in
        I)      INP="$OPTARG";;
        L)      CHR="$OPTARG";;
        R)      REF="$OPTARG";;
        D)      DBSNP="$OPTARG";;
        E)      ExonFile="$OPTARG";;
	h)	echo $USAGE
		exit 1;;
	[?])	echo $USAGE
		exit 1;;
        esac
done

if [[ $INP == "" || $REF == "" || $DBSNP == "" ]]
then
	echo $USAGE
	exit 1
fi


echo
date
echo "Calling ./gatk_realign.scr step by step..."

if [[ $CHR == "" ]]
then
	qsub -l mem=8G,time=16:: -o realignment.X.output -e realignment.X.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L "X"
	qsub -l mem=8G,time=16:: -o realignment.Y.output -e realignment.Y.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L "Y"
	qsub -l mem=8G,time=16:: -t 1-22 -sync y ./gatk_realign.scr -I $INP -R $REF -D $DBSNP
#	qsub -l mem=8G,time=16:: -o realignment.$SGE_TASK_ID.output -e realignment.$SGE_TASK_ID.output -t 1-22 -sync y ./gatk_realign.scr -I $INP -R $REF -D $DBSNP

	mv realignment.1.output realignment.output
	mv $INP.1.fixed.bam $INP.fixed.bam
	mv $INP.1.fixed.bam.bai $INP.fixed.bam.bai
	for i in 1 2 3 4 5 6 7 8 9 X Y 10 11 12 13 14 15 16 17 18 19 20 21 22
	do
        	cat realignment.$i.output >> realignment.output
        	rm realignment.$i.output
        	cat $INP.$i.fixed.bam >> $INP.fixed.bam
        	rm $INP.$i.fixed.bam
        	cat $INP.$i.fixed.bam.bai >> $INP.fixed.bam.bai
        	rm $INP.$i.fixed.bam.bai
	done
else
	qsub -l mem=8G,time=16:: -sync y ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L $CHR
fi
if [[ $? != 0 || `grep "$ERRORMESSAGE" realignment.output pipeline.output` != "" ]]
then
	echo "Realignment FAILED"
	exit 1
fi
echo "Execution of ./gatk_realign.scr complete"
echo



date
echo "Calling ./gatk_calibrate.scr ..."
qsub -l mem=8G,time=8:: -sync y ./gatk_calibrate.scr -I $INP -L $CHR -R $REF -D $DBSNP 2>&1 > calibration.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" calibration.output pipeline.output` != "" ]]
then
	echo "Calibration FAILED"
	exit 1
fi
echo "Execution of ./gatk_calibrate.scr complete"
echo



date
echo "Calling ./gatk_recalibrate.scr ..."
qsub -l mem=8G,time=8:: -sync y ./gatk_recalibrate.scr -I $INP -L $CHR -R $REF 2>&1 > recalibration.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" recalibration.output pipeline.output` != "" ]]
then
	echo " Recalibration FAILED"
	exit 1
fi
echo "Execution of ./gatk_recalibrate.scr complete"
echo



date
echo "Calling ./gatk_depthofcoverage.scr ..."
if [[ $ExonFile == "" ]]
then
	qsub -l mem=8G,time=8:: ./gatk_depthofcoverage.scr -I $INP -L $CHR -R $REF 2>&1 > depthofcoverage.output
else
	qsub -l mem=8G,time=8:: ./gatk_depthofcoverage.scr -I $INP -E "$ExonFile" -R $REF 2>&1 > depthofcoverage.output
fi
if [[ $? != 0 || `grep "$ERRORMESSAGE" depthofcoverage.output pipeline.output` != "" ]]
then
        echo "Depth of coverage FAILED"
        exit 1
fi
echo "Execution of ./gatk_depthofcoverage.scr complete"
echo



date
echo "Calling ./gatk_indelcall.scr ..."
qsub -l mem=8G,time=8:: -sync y ./gatk_indelcall.scr -I $INP -L $CHR -R $REF 2>&1 > indelcalling.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" indelcalling.output pipeline.output` != "" ]]
then
	echo "Indel calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_indelcall.scr complete"
echo



date
echo "Calling ./gatk_snpcall.scr ..."
qsub -l mem=8G,time=8:: -sync y ./gatk_snpcall.scr -I $INP -L $CHR -R $REF -D $DBSNP 2>&1 > snpcalling.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" snpcalling.output pipeline.output` != "" ]]
then
	echo "SNP calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_snpcall.scr complete"
echo



date
echo "Calling ./gatk_varianteval.scr ..."
qsub -l mem=1G,time=:30: -sync y ./gatk_varianteval.scr -I $INP -L $CHR -R $REF -D $DBSNP 2>&1 > varianteval_prefilter.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" varianteval_prefilter.output pipeline.output` != "" ]]
then
	echo "Prefilter variant eval calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_varianteval.scr complete"
echo



date
echo "Calling ./gatk_filter.scr ..."
qsub -l mem=1G,time=:30: -sync y ./gatk_filter.scr -I $INP -L $CHR -R $REF 2>&1 > variantfiltering.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" variantfiltering.output pipeline.output` != "" ]]
then
	echo "Variant filtering FAILED"
	exit 1
fi
echo "Execution of ./gatk_filter.scr complete"
echo



date
echo "Calling ./gatk_varianteval.scr ..."
qsub -l mem=1G,time=:30: -sync y ./gatk_varianteval.scr -I $INP -L $CHR -R $REF -D $DBSNP 2>&1 > varianteval_postfilter.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" varianteval_postfilter.output pipeline.output` != "" ]]
then
	echo "Postfilter variant eval calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_varianteval.scr complete"
echo
date

