#!/bin/bash
#$ -cwd
# Findmem

# Realign reads locally around indels
# Run with all samples from a single chromosome together in $INP.list

HEAP=4000
TEMP="/ifs/scratch/c2b2/af_lab/saec/temp/"
INP=""
CHR=""
REF=""
DBSNP=""
RMODE="0"  # default is not in batch mode
maxsize=130000000
SIZE=0 

USAGE="Usage: $0 -I <Input bam file> -R <Reference fasta> -D <DBSNP file> [-L \"#:#-#\"]"
ERRORMESSAGE="#### ERROR"
ERRORMESSAGE1="The following error has occurred"

while getopts I:L:R:D:g:m:s:h o
do      case "$o" in
        I)      INP="$OPTARG";;
        L)      CHR="$OPTARG";;
        R)      REF="$OPTARG";;
        D)      DBSNP="$OPTARG";;
        g)      GLOBAL="$OPTARG";;  # global config
        m)      RMODE="$OPTARG";;  # running mode; none-zero means running in batch mode
        s)      SIZE="$OPTARG";;
        h)      echo $USAGE
                exit 1;;
        esac
done


if [[ $INP == ""  || $REF == "" || $DBSNP == "" || $GLOBAL == "" ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL

JAVA="java -Xmx${HEAP}m -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}

# FIXMATE="$JAVA -jar /ifs/data/c2b2/ip_lab/shares/SOFTWARE/picard-tools-1.29/FixMateInformation.jar"
FIXMATE="$JAVA -jar ${PICARD}/FixMateInformation.jar" 

# /ifs/home/c2b2/af_lab/saec/Software/picard-tools-1.33/FixMateInformation.jar"
# SAMTOOLS="/ifs/home/c2b2/ip_lab/yshen/usr/bin/samtools"

OUTDIR=$INP"_pipe"
if [ ! -d $OUTDIR ]; then
  mkdir $OUTDIR
fi


# split a chromosome into two if the size is larger than 170M
if [[ $RMODE == "0" ]]  # running in non-batch mode 
    then
    if [[ $CHR == "" ]]
	then
        CHR="${SGE_TASK_ID}"
    fi
    if [[ $CHR == "23" ]]
	then
	CHR="X"
    fi
    if [[ $CHR == "24"  ]]
	then
	CHR="Y"
    fi
    if [[ $REFTYPE == "hg" ]]  # hg18/19 -> chr1, chr2 etc;  build36/37 -> 1, 2 etc                                                                                                                       
	then
        CHR="chr${CHR}"
    fi
    
    if [ ! -e $INP.header ]
    then # 
      ${SAMTOOLS} view -H $INP -o $INP.header 
    fi
    SIZE=`grep -w SN:$CHR $INP.header | cut -f3 | sed 's/LN\://'`
    echo "size of the chromosome $CHR:  $SIZE"
    if [[ $SIZE -gt $maxsize ]]
	then
#
#	newchr1="$CHR:1-$halfs"
#	newchr2="$CHR:$halfs-$chrsize"


	qsub -l mem=5G,time=24:: -t 1-2 -sync y -o $OUTDIR/$CHR.realign.stdout -e $OUTDIR/$CHR.realign.stderr ${BPATH}/gatk_realign.scr  -g $GLOBAL -I $INP -R $REF -D $DBSNP -m 1 -L $CHR -s $SIZE


	
	exit 0 
    fi
else  # run in batch model
    halfs=$(($SIZE/2 + 1))
    	# echo "chr: $CHR ; size: $SIZE"
    if [[  ${SGE_TASK_ID} == "1" ]] # first half
	then
	CHR="$CHR:1-$halfs"

    else
	CHR="$CHR:$halfs-$SIZE"
	# echo "chr: $CHR"
    fi
fi

echo $CHR 
function realigner {



}

$GATK \
 -L $CHR \
 -T RealignerTargetCreator \
 -I $INP \
 -R $REF \
 -D $DBSNP \
 -o $OUTDIR/$CHR.forRealigner.intervals 2>&1 >> $OUTDIR/$CHR.realign.out

$GATK \
 -L $CHR \
 -I $INP \
 -R $REF \
 -D $DBSNP \
 -T IndelRealigner \
 -compress 0 \
 -targetIntervals $OUTDIR/$CHR.forRealigner.intervals \
 -o $OUTDIR/$CHR.cleaned.bam 2>&1 >> $OUTDIR/$CHR.realign.out

# Need to implement for known-only/lane level cleaning?


$FIXMATE INPUT=$OUTDIR/$CHR.cleaned.bam OUTPUT=$OUTDIR/$CHR.fixed.bam SO=coordinate VALIDATION_STRINGENCY=SILENT 


echo "$0: realign complete"

rm -f $OUTDIR/$CHR.cleaned.bam
rm -rf $OUTDIR/$CHR.realign.*
