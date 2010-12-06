#!/bin/sh
#$ -cwd

GLOBAL="global_config.sh"

if [[ -e $GLOBAL ]]
then
        . $GLOBAL
else
	echo "Global config file not found. Exiting."
	exit 1
fi

USAGE="Usage: run_pipeline.sh -I <Input bam file> -R <Reference fasta> -D <DBSNP file> -P <Platform> [[-L #:#-#] | [-E <ExonFile>]]"

while getopts I:L:R:D:P:E:h o
do      case "$o" in
        I)      INP="$OPTARG";;
        L)      CHR="$OPTARG";;
        R)      REF="$OPTARG";;
        D)      DBSNP="$OPTARG";;
        P)      Platform="$OPTARG";;
        E)      ExonFile="$OPTARG";;
	h)	echo $USAGE
		exit 1;;
        esac
done

if [[ $INP == "" || $REF == "" || $DBSNP == "" || $Platform == "" ]]
then
	echo $USAGE
	exit 1
fi

SUFFIX=".fixed.bam"
SUFFIX_BAI=".fixed.bam.bai"
DATAPATH=`dirname "$INP"`
CONTIG_ORDER="`$SAMTOOLS idxstats $INP | grep -m 24 [0-9] | awk '{print $1}' | tr '\n' ' '`"
CONTIG_INP_ORDER="`$SAMTOOLS idxstats $INP | grep -m 24 [0-9] | awk -v var1="$INP" -v var2="$SUFFIX" '{print var1 $1 var2}' | tr '\n' ' '`"
CONTIG_INP_BAI_ORDER="`$SAMTOOLS idxstats $INP | grep -m 24 [0-9] | awk -v var1="$INP" -v var2="$SUFFIX_BAI" '{print var1 $1 var2}' | tr '\n' ' '`"
RESULT=0

if [[ $CHR == "" && $ExonFile == "" ]]
then
	echo $EXONUSAGE
	exit 1
fi

realignment_count=0
while [[ $RESULT == 0 ]]
do
	date
	realignment_count=`expr $realignment_count + 1`
	rm -f $DATAPATH/realignment*.output
	echo Calling ./gatk_realign.scr step by step... Attempt "$realignment_count"
	RESULT=1
	if [[ $CHR == "" ]]
	then
		RETX=`qsub -l mem=5G,time=2:: -e $DATAPATH/realignment.output -o $DATAPATH/realignment.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L "X"`
		sleep 10
		RETY=`qsub -l mem=5G,time=2:: -e $DATAPATH/realignment.output -o $DATAPATH/realignment.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L "Y"`
		sleep 10
		RET=`qsub -l mem=6G,time=16:: -t 1-22 -e $DATAPATH/realignment.output -o $DATAPATH/realignment.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP | awk '{print $3}'`
	else
		RET=`qsub -l mem=8G,time=16:: -o $DATAPATH/realignment.output -e $DATAPATH/realignment.output ./gatk_realign.scr -I $INP -R $REF -D $DBSNP -L $CHR | awk '{print $3}'`
	fi

	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running realigner...
		sleep 100
		if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/realignment*.output "$DATAPATH"/pipeline.output` != "" ]]
		then
			echo "Realignment FAILED"
			qdel $RET $RETX $RETY
			cp $DATAPATH/pipeline.output $DATAPATH/pipeline.output.$realignment_count
			cp $DATAPATH/realignment.output $DATAPATH/realignment.output.$realignment_count
			> $DATAPATH/pipeline.output
			> $DATAPATH/realignment.output
			RESULT=0
			break
		fi
	done
done
echo "Execution of ./gatk_realign.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/merge.output
	echo "Calling ./gatk_merge.scr ..."
	RESULT=1
	if [[ $CHR == "" ]]
	then
		RET=`qsub -l mem=2G,time=4:: -o $DATAPATH/merge.output -e $DATAPATH/merge.output ./gatk_merge.scr -I $INP.fixed.bam -R $REF -C "$CONTIG_ORDER" -O "$CONTIG_INP_ORDER" -B "$CONTIG_INP_BAI_ORDER" | awk '{print $3}'`
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		while [[ -z $STATUS ]]
		do
			STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
			echo Running merge operation...
			sleep 10
		done
		if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/merge.output "$DATAPATH"/pipeline.output` != "" ]]
		then
        		echo "Merging of realigned bam files FAILED"
			RESULT=0
		fi
	fi
done
echo "Execution of ./gatk_merge.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/calibrate.output
	echo "Calling ./gatk_calibrate.scr ..."
	RESULT=1
	RET=`qsub -l mem=8G,time=8:: -o $DATAPATH/calibrate.output -e $DATAPATH/calibrate.output ./gatk_calibrate.scr -I $INP.fixed.bam -R $REF -D $DBSNP -L $CHR | awk '{print $3}'`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running Calibration...
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/calibrate.output "$DATAPATH"/pipeline.output` != "" ]]
	then
		echo "Calibration FAILED"
		RESULT=0
	fi
done
echo "Execution of ./gatk_calibrate.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/recalibrate.output
	echo "Calling ./gatk_recalibrate.scr ..."
	RESULT=1
	RET=`qsub -l mem=8G,time=8:: -o $DATAPATH/recalibrate.output -e $DATAPATH/recalibrate.output ./gatk_recalibrate.scr -I $INP.fixed.bam -R $REF -L $CHR | awk '{print $3}'`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running Recalibration...
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/recalibrate.output "$DATAPATH"/pipeline.output` != "" ]]
	then
		echo " Recalibration FAILED"
		RESULT=0
	fi
done
echo "Execution of ./gatk_recalibrate.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/depthofcoverage.output
	echo "Calling ./gatk_depthofcoverage.scr ..."
	RESULT=1
	if [[ $ExonFile == "" ]]
	then
		qsub -l mem=6G,time=8:: -o $DATAPATH/depthofcoverage.output -e $DATAPATH/depthofcoverage.output ./gatk_depthofcoverage.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -L $CHR
	else
		qsub -l mem=6G,time=8:: -o $DATAPATH/depthofcoverage.output -e $DATAPATH/depthofcoverage.output ./gatk_depthofcoverage.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -E $ExonFile
	fi

	while [[ ! -s $DATAPATH/depthofcoverage.output ]]
	do
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/depthofcoverage.output "$DATAPATH"/pipeline.output` != "" ]]
	then
	        echo "Depth of coverage FAILED"
		echo Running Depth Coverage...
		RESULT=0
	fi
done
echo "Execution of ./gatk_depthofcoverage.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/indelcalling.output
	echo "Calling ./gatk_indelcall.scr ..."
	RESULT=1
	RET=`qsub -l mem=6G,time=8:: -o $DATAPATH/indelcalling.output -e $DATAPATH/indelcalling.output ./gatk_indelcall.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -L $CHR | awk '{print $3}'`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running Indel Calling...
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/indelcalling.output "$DATAPATH"/pipeline.output` != "" ]]
	then
		echo "Indel calling FAILED"
		RESULT=0
	fi
done
echo "Execution of ./gatk_indelcall.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/snpcalling.output
	echo "Calling ./gatk_snpcall.scr ..."
	RESULT=1
	RET=`qsub -l mem=6G,time=8:: -o $DATAPATH/snpcalling.output -e $DATAPATH/snpcalling.output ./gatk_snpcall.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -D $DBSNP -P $Platform -L $CHR | awk '{print $3}'`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running SNP Calling...
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/snpcalling.output "$DATAPATH"/pipeline.output` != "" ]]
	then
		echo "SNP calling FAILED"
		RESULT=0
	fi
done
echo "Execution of ./gatk_snpcall.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/varianteval_prefiltering.output
	echo "Calling ./gatk_varianteval.scr ..."
	RESULT=1
	RET=`qsub -l mem=5G,time=1:: -o $DATAPATH/varianteval_prefiltering.output -e $DATAPATH/varianteval_prefiltering.output ./gatk_varianteval.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -D $DBSNP | awk '{print $3}'`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running Prefilter Variant Evaluation...
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/varianteval_prefiltering.output "$DATAPATH"/pipeline.output` != "" ]]
	then
		echo "Prefilter variant eval calling FAILED"
		RESULT=0
	fi
done
echo "Execution of ./gatk_varianteval.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/variantfiltering.output
	echo "Calling ./gatk_filter.scr ..."
	RESULT=1
	RET=`qsub -l mem=5G,time=1:: -o $DATAPATH/variantfiltering.output -e $DATAPATH/variantfiltering.output ./gatk_filter.scr -I $INP.fixed.bam.recalibrated.bam -R $REF | awk '{print $3}'`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running Variant Filtration...
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/variantfiltering.output "$DATAPATH"/pipeline.output` != "" ]]
	then
		echo "Variant filtering FAILED"
		RESULT=0
	fi
done
echo "Execution of ./gatk_filter.scr complete"
RESULT=0
echo



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/varianteval_postfiltering.output
	echo "Calling ./gatk_varianteval.scr ..."
	RESULT=1
	RET=`qsub -l mem=5G,time=1:: -o $DATAPATH/varianteval_postfiltering.output -e $DATAPATH/varianteval_postfiltering.output ./gatk_varianteval.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -D $DBSNP | awk '{print $3}'`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running Postfilter Variant Evaluation
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/varianteval_postfiltering.output "$DATAPATH"/pipeline.output` != "" ]]
	then
		echo "Postfilter variant eval calling FAILED"
		RESULT=0
	fi
done
echo "Execution of ./gatk_varianteval.scr complete"
RESULT=0
echo
date



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/compute_ratios.output
	echo "Calling ./gatk_compute_ratios.scr ..."
	RESULT=1
	RET=`qsub -l mem=5G,time=4:: -o $DATAPATH/compute_ratios.output -e $DATAPATH/compute_ratios.output ./gatk_compute_ratios.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -D $DBSNP | awk '{print $3}'`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running Ratio Compute job
		sleep 10
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" "$DATAPATH"/compute_ratios.output "$DATAPATH"/pipeline.output` != "" ]]
	then
		echo "Compute Ratio job FAILED"
		RESULT=0
	fi
done
echo "Execution of ./gatk_compute_ratios.scr complete"
RESULT=0
echo
date



while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/count_snv.output
	echo "Calling ./count_all_SNV.sh ..."
	RESULT=1
	RET=`qsub -l mem=5G,time=4:: -o $DATAPATH/count_snv.output -e $DATAPATH/count_snv.output ./count_all_SNV.sh -I $INP.fixed.bam.recalibrated.bam -V $INP.fixed.bam.recalibrated.bam.snps.filtered.vcf -E $ExonFile`
	STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
	while [[ -z $STATUS ]]
	do
		STATUS=`qstat -j "$RET" 2>&1 | grep "$JOB_STAT"`
		echo Running SNV Count job
		sleep 30
	done
	if [[ $? != 0 || `grep "$ERRORMESSAGE" $DATAPATH/count_snv.output $DATAPATH/pipeline.output` != "" ]]
	then
		echo Running SNV Count job FAILED
		RESULT=0
	fi
done
echo "Execution of ./count_all_SNV.sh complete"
RESULT=0
echo
date

chmod g=rw $DATAPATH/*

