#!/bin/sh
#$ -cwd

unset $CHR 

INP=""
CHR=""
REF=""
DBSNP=""
Platform=""
refseq=""
ExonFile=""

USAGE="Usage: $0 -I <Input bam file> -R <Reference fasta> -D <DBSNP file> -S <refseq> -P <Platform> [[-L #:#-#] | [-E <ExonFile>]]"
ERRORMESSAGE="The following error has occurred"

while getopts I:L:R:D:P:S:E:h o
do      case "$o" in
        I)      INP="$OPTARG";;
        L)      CHR="$OPTARG";;
        R)      REF="$OPTARG";;
        D)      DBSNP="$OPTARG";;
        P)      Platform="$OPTARG";;
        S)      refseq="$OPTARG";;
        E)      ExonFile="$OPTARG";;
	h)	echo $USAGE
		exit 1;;
	[?])	echo $USAGE
		exit 1;;
        esac
done

if [[ $INP == "" || $REF == "" || $DBSNP == "" || $refseq == "" || $Platform == "" ]]
then
	echo $USAGE
	exit 1
fi
if [[ $CHR != "" && $ExonFile != "" ]]
then
	echo $USAGE
	exit 1
fi


echo
echo "Calling ./gatk_depthofcoverage.scr ..."
if [[ $ExonFile == "" ]]
then
	./gatk_depthofcoverage.scr -I $INP -L $CHR -R $REF 2>&1 > depthofcoverage.output
else
	cat $ExonFile | awk '{split($1,a,"chr"); print a[2] ":" $2 "-" $3}' > "${ExonFile}.list"
	./gatk_depthofcoverage.scr -I $INP -L "$ExonFile.list" -R $REF 2>&1 > depthofcoverage.output
fi
if [[ $? != 0 || `grep "$ERRORMESSAGE" depthofcoverage.output` != "" ]]
then
        echo "Depth of coverage FAILED"
        exit 1
fi
echo "Execution of ./gatk_depthofcoverage.scr complete"
echo



echo "Calling ./gatk_calibrate.scr ..."
./gatk_calibrate.scr -I $INP -L $CHR -R $REF -D $DBSNP -P $Platform 2>&1 > calibration.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" calibration.output` != "" ]]
then
	echo "Calibration FAILED"
	exit 1
fi
echo "Execution of ./gatk_calibrate.scr complete"
echo



echo "Calling ./gatk_recalibrate.scr ..."
./gatk_recalibrate.scr -I $INP -L $CHR -R $REF -P $Platform 2>&1 > recalibration.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" recalibration.output` != "" ]]
then
	echo " Recalibration FAILED"
	exit 1
fi
echo "Execution of ./gatk_recalibrate.scr complete"
echo



echo "Calling ./gatk_realign.scr ..."
./gatk_realign.scr -I $INP -L $CHR -R $REF -D $DBSNP 2>&1 > realignment.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" realignment.output` != "" ]]
then
	echo "Realignment FAILED"
	exit 1
fi
echo "Execution of ./gatk_realign.scr complete"
echo



echo "Calling ./gatk_indelcall.scr ..."
./gatk_indelcall.scr -I $INP -L $CHR -R $REF -S $refseq 2>&1 > indelcalling.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" indelcalling.output` != "" ]]
then
	echo "Indel calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_indelcall.scr complete"
echo




echo "Calling ./gatk_snpcall.scr ..."
./gatk_snpcall.scr -I $INP -L $CHR -R $REF -D $DBSNP -P $Platform 2>&1 > snpcalling.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" snpcalling.output` != "" ]]
then
	echo "SNP calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_snpcall.scr complete"
echo


echo "Calling ./gatk_varianteval.scr ..."
./gatk_varianteval.scr -I $INP -L $CHR -R $REF -D $DBSNP 2>&1 > varianteval_prefilter.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" varianteval_prefilter.output` != "" ]]
then
	echo "Prefilter variant eval calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_varianteval.scr complete"
echo



echo "Calling ./gatk_filter.scr ..."
./gatk_filter.scr -I $INP -L $CHR -R $REF 2>&1 > snpcalling.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" variantfiltering.output` != "" ]]
then
	echo "Variant filtering FAILED"
	exit 1
fi
echo "Execution of ./gatk_filter.scr complete"
echo




echo "Calling ./gatk_varianteval.scr ..."
./gatk_varianteval.scr -I $INP -L $CHR -R $REF -D $DBSNP 2>&1 > varianteval_postfilter.output
if [[ $? != 0 || `grep "$ERRORMESSAGE" varianteval_postfilter.output` != "" ]]
then
	echo "Postfilter variant eval calling FAILED"
	exit 1
fi
echo "Execution of ./gatk_varianteval.scr complete"
echo

