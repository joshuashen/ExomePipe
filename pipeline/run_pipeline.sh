#!/bin/sh
#$ -cwd


# The main script that invokes all steps part of the pipeline. The order of the steps is fixed. If any intermediate steps fails, the next step will not be executed.
# This is true for all but the step to compute the depth of coverage. The depth of coverage computation is kicked off in parallel to the other steps in the pipeline
# since no other step is dependent on it.
# An automatic resubmission of jobs logic has been built. It may not be fail proof, but it works most of the time thus saving much manual intervention.
# The date command is invoked before and after each step so as to find out the time taken for executing that step.

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


# The CONTIG_* variables give information on the contig order of the chromosomes for the given bam file under consideration.
# Based on this order, they contain a generated string which is used to perform a specific action.
# This is useful in the post realignment bam file merge step.
SUFFIX=".fixed.bam"
SUFFIX_BAI=".fixed.bam.bai"
DATAPATH=`dirname "$INP"`
CONTIG_ORDER="`$SAMTOOLS idxstats $INP | grep -m 24 [0-9] | awk '{print $1}' | tr '\n' ' '`"
CONTIG_INP_ORDER="`$SAMTOOLS idxstats $INP | grep -m 24 [0-9] | awk -v var1="$INP" -v var2="$SUFFIX" '{print var1 $1 var2}' | tr '\n' ' '`"
CONTIG_INP_BAI_ORDER="`$SAMTOOLS idxstats $INP | grep -m 24 [0-9] | awk -v var1="$INP" -v var2="$SUFFIX_BAI" '{print var1 $1 var2}' | tr '\n' ' '`"
RESULT=0


# Either or both of the chromosome # or the file containing exon region details needs to be provided.
# Normally, one of the two is to be provided and not both.
if [[ $CHR == "" && $ExonFile == "" ]]
then
	echo $EXONUSAGE
	exit 1
fi


# This steps performs realignment of the bam file. It is the most time consuming step and so it is split into smaller jobs.
# The while loop aids in automatic re-running of failed jobs. A failed job is determined based on the status returned by the qstat command as described below.
# The realignment_count variable keeps track of the number of attempts. The logs from each attempt are saved based on this number.
realignment_count=0
while [[ $RESULT == 0 ]]
do
	date
	realignment_count=`expr $realignment_count + 1`
	rm -f $DATAPATH/realignment*.output
	echo Calling ./gatk_realign.scr step by step... Attempt "$realignment_count"
	RESULT=1

	# If chromosome # is not specified it is assumed that realignment needs to be run for the entire chromosome. This step is then split up to run for the
	# X, Y and chromosomes 1-22. If chromosome # was specified, realignment is run for that chromosome alone in the else condition.
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

	# The status of the job is checked periodically until the job has ended. When the job has ended, STATUS will be non-empty since grep will return a value.
	# There are two ways we check for an error:
	#   1. An error has occured in the environment and $? is non-zero.
	#   2. An error message is printed in either of the realignment or the run_pipeline log file.
	# If either of the two is found, all remaining split-jobs are killed, their logs saved and we loop again.
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


# Post realignment, the individual realigned bam files are merged. This is done only if a chromosome # was not specified or else there were no split-jobs.
# It is necessary to merge the bam files before the next step of calibration since calibration gives more accurate results when performed on the entire genome.
# The looping logic is similar to the realignment step above.
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


# This step performs calibration of the merged bam file from the previous step.
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


# This step performs recalibration of the bam file from the previous step.
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


# This step computes the depth of coverage obtained. It is important to note that this step is performed on the recalibrated file and not on the original bam file.
# Also while there is an auto resubmit logic built for this step as well, failure of this step does not block the remaining pipeline from running. In fact all we
# wait for is the step to begin running before we proceed to the next step.
while [[ $RESULT == 0 ]]
do
	date
	rm -f $DATAPATH/depthofcoverage.output
	echo "Calling ./gatk_depthofcoverage.scr ..."
	RESULT=1

	# If both chromosome # and the Exon list are specified, priority is given to the Exon list. Hence depth of coverage will be computed for the entire set of
	# exon regions. Normally we would expect one of the two parameters alone to be specified.
	if [[ $ExonFile == "" ]]
	then
		qsub -l mem=6G,time=8:: -o $DATAPATH/depthofcoverage.output -e $DATAPATH/depthofcoverage.output ./gatk_depthofcoverage.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -L $CHR
	else
		qsub -l mem=6G,time=8:: -o $DATAPATH/depthofcoverage.output -e $DATAPATH/depthofcoverage.output ./gatk_depthofcoverage.scr -I $INP.fixed.bam.recalibrated.bam -R $REF -E $ExonFile
	fi

	# We wait until the qsub command has actually kicked in. When this happens, the output log file will be created.
	while [[ ! -s $DATAPATH/depthofcoverage.output ]]
	do
		sleep 10
	done

	# For this step we check for failure only towards the beginning of execution. It is quite complicated to check for failure in a loop until completion.
	# Hence there may be a little more manual intervention involved in this step. However, a good number of times if the step was not invoked correctly etc.,
	# the step fails in the very beginning and is caught. While this is reported in the log file, we still permit execution for the rest of the script.
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


# This step performs indel calling on the recalibrated bam file.
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


# This step performs SNP calling on the recalibrated bam file.
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


# This step performs variant evaluation on the recalibrated bam file.
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


# This step performs variant filtering on the recalibrated bam file.
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


# This step performs variant evaluation on the recalibrated bam file.
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


# This step computes the transition/transversion ratio.
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


# This step performs am SNV count for the filtered vcfs.
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

# Provide read/write access to all generated files.
chmod g=rw $DATAPATH/*

