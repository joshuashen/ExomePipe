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
	qsub -l mem=8G,time=6:: -o realignment.X.output -e realignment.X.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L "X"
	qsub -l mem=8G,time=6:: -o realignment.Y.output -e realignment.Y.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L "Y"
	qsub -l mem=8G,time=16:: -t 1-22 -sync y -o realignment.output -e realignment.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP
else
	qsub -l mem=8G,time=16:: -sync y -o realignment.output -e realignment.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L $CHR
fi
if [[ $? != 0 || `grep "$ERRORMESSAGE" realignment.*.output pipeline.output` != "" ]]
then
	echo "Realignment FAILED"
	exit 1
fi
echo "Execution of ./gatk_realign.scr complete"
echo



date
echo "Calling ./gatk_calibrate.scr ..."
if [[ $CHR == "" ]]
then
	qsub -l mem=8G,time=6:: -o calibrate.X.output -e calibrate.X.output ./gatk_calibrate.scr -I $INP -R $REF -D $DBSNP -L "X"
	qsub -l mem=8G,time=6:: -o calibrate.Y.output -e calibrate.Y.output ./gatk_calibrate.scr -I $INP -R $REF -D $DBSNP -L "Y"
	qsub -l mem=8G,time=16:: -t 1-22 -sync y -o calibrate.output -e calibrate.output ./gatk_calibrate.scr -I $INP -R $REF -D $DBSNP
else
	qsub -l mem=8G,time=16:: -sync y -o calibrate.output -e calibrate.output ./gatk_calibrate.scr -I $INP -R $REF -D $DBSNP -L $CHR
fi
if [[ $? != 0 || `grep "$ERRORMESSAGE" calibration.*.output pipeline.output` != "" ]]
then
	echo "Calibration FAILED"
	exit 1
fi
echo "Execution of ./gatk_calibrate.scr complete"
echo



date
echo "Calling ./gatk_recalibrate.scr ..."
if [[ $CHR == "" ]]
then
	qsub -l mem=8G,time=6:: -o recalibrate.X.output -e recalibrate.X.output ./gatk_recalibrate.scr -I $INP -R $REF -L "X"
	qsub -l mem=8G,time=6:: -o recalibrate.Y.output -e recalibrate.Y.output ./gatk_recalibrate.scr -I $INP -R $REF -L "Y"
	qsub -l mem=8G,time=16:: -t 1-22 -sync y -o recalibrate.output -e recalibrate.output ./gatk_recalibrate.scr -I $INP -R $REF
else
	qsub -l mem=8G,time=16:: -sync y -o recalibrate.output -e recalibrate.output ./gatk_recalibrate.scr -I $INP -R $REF -L $CHR
fi
if [[ $? != 0 || `grep "$ERRORMESSAGE" recalibration.*.output pipeline.output` != "" ]]
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
	qsub -l mem=8G,time=8:: -o depthofcoverage.output -e depthofcoverage.output ./gatk_depthofcoverage.scr -I $INP -R $REF -L $CHR
else
	qsub -l mem=8G,time=8:: -o depthofcoverage.X.output -e depthofcoverage.X.output ./gatk_depthofcoverage.scr -I $INP -R $REF -E "$ExonFile" -L "X"
	qsub -l mem=8G,time=8:: -o depthofcoverage.Y.output -e depthofcoverage.Y.output ./gatk_depthofcoverage.scr -I $INP -R $REF -E "$ExonFile" -L "Y"
	qsub -l mem=8G,time=8:: -t 1-22 -o depthofcoverage.output -e depthofcoverage.output ./gatk_depthofcoverage.scr -I $INP -R $REF -E "$ExonFile"
fi
if [[ $? != 0 || `grep "$ERRORMESSAGE" depthofcoverage.*.output pipeline.output` != "" ]]
then
        echo "Depth of coverage FAILED"
        exit 1
fi
echo "Execution of ./gatk_depthofcoverage.scr complete"
echo



date
echo "Calling ./gatk_indelcall.scr ..."
if [[ $CHR == "" ]]
then
	qsub -l mem=8G,time=6:: -o indelcalling.X.output -e indelcalling.X.output ./gatk_indelcall.scr -I $INP -R $REF -L "X"
	qsub -l mem=8G,time=6:: -o indelcalling.Y.output -e indelcalling.Y.output ./gatk_indelcall.scr -I $INP -R $REF -L "Y"
	qsub -l mem=8G,time=16:: -t 1-22 -sync y -o indelcalling.output -e indelcalling.output ./gatk_indelcall.scr -I $INP -R $REF
else
	qsub -l mem=8G,time=16:: -sync y -o indelcalling.output -e indelcalling.output ./gatk_indelcall.scr -I $INP -R $REF -L $CHR
fi
mv $INP.1.indels.raw.bed $INP.indels.raw.bed
mv $INP.1.detailed.output.bed $INP.detailed.output.bed
mv $INP.1.output.vcf $INP.output.vcf
for i in 2 3 4 5 6 7 8 9 X Y 10 11 12 13 14 15 16 17 18 19 20 21 22
do
	cat $INP.$i.indels.raw.bed >> $INP.indels.raw.bed
	rm $INP.$i.indels.raw.bed
	cat $INP.$i.detailed.output.bed >> $INP.detailed.output.bed
	rm $INP.$i.detailed.output.bed
	cat $INP.$i.output.vcf >> $INP.output.vcf
	rm $INP.$i.output.vcf
done
if [[ $? != 0 || `grep "$ERRORMESSAGE" indelcalling.*.output pipeline.output` != "" ]]
then
	echo "Indel calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_indelcall.scr complete"
echo



date
echo "Calling ./gatk_snpcall.scr ..."
if [[ $CHR == "" ]]
then
	qsub -l mem=8G,time=6:: -o snpcalling.X.output -e snpcalling.X.output ./gatk_snpcall.scr -I $INP -R $REF -D $DBSNP -L "X"
	qsub -l mem=8G,time=6:: -o snpcalling.Y.output -e snpcalling.Y.output ./gatk_snpcall.scr -I $INP -R $REF -D $DBSNP -L "Y"
	qsub -l mem=8G,time=16:: -t 1-22 -sync y -o snpcalling.output -e snpcalling.output ./gatk_snpcall.scr -I $INP -R $REF -D $DBSNP
else
	qsub -l mem=8G,time=16:: -sync y -o snpcalling.output -e snpcalling.output ./gatk_snpcall.scr -I $INP -R $REF -D $DBSNP -L $CHR
fi
if [[ $? != 0 || `grep "$ERRORMESSAGE" snpcalling.*.output pipeline.output` != "" ]]
then
	echo "SNP calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_snpcall.scr complete"
echo



date
echo "Calling ./gatk_merge.scr ..."
qsub -l mem=2G,time=1:: -sync y -o merge.output -e merge.output ./gatk_merge.scr -I $INP -R $REF
if [[ $? != 0 || `grep "$ERRORMESSAGE" merge.output pipeline.output` != "" ]]
then
        echo "Merging of vcf files FAILED"
        exit 1
fi
echo "Execution of ./gatk_merge.scr complete"
echo



date
echo "Calling ./gatk_varianteval.scr ..."
qsub -l mem=5G,time=1:: -sync y -o varianteval_prefiltering.output -e varianteval_prefiltering.output ./gatk_varianteval.scr -I $INP -R $REF -D $DBSNP
if [[ $? != 0 || `grep "$ERRORMESSAGE" varianteval_prefiltering.output pipeline.output` != "" ]]
then
	echo "Prefilter variant eval calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_varianteval.scr complete"
echo



date
echo "Calling ./gatk_filter.scr ..."
qsub -l mem=5G,time=1:: -sync y -o variantfiltering.output -e variantfiltering.output ./gatk_filter.scr -I $INP -R $REF
if [[ $? != 0 || `grep "$ERRORMESSAGE" variantfiltering.output pipeline.output` != "" ]]
then
	echo "Variant filtering FAILED"
	exit 1
fi
echo "Execution of ./gatk_filter.scr complete"
echo



date
echo "Calling ./gatk_varianteval.scr ..."
qsub -l mem=5G,time=1:: -sync y -o varianteval_postfiltering.output -e varianteval_postfiltering.output ./gatk_varianteval.scr -I $INP -R $REF -D $DBSNP
if [[ $? != 0 || `grep "$ERRORMESSAGE" varianteval_postfiltering.output pipeline.output` != "" ]]
then
	echo "Postfilter variant eval calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_varianteval.scr complete"
echo
date

