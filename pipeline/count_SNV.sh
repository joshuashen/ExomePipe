#!/bin/sh
#$ -cwd

CHECK1="HaplotypeScore"
CHECK2="HRun"

USAGE="Usage: SNV_Count.sh -V <Variants file> -E <ExonFile> -L <CHR>"

while getopts V:E:L:h o
do      case "$o" in
        V)      VCF_FILE="$OPTARG";;
        E)      EXON_FILE="$OPTARG";;
        L)      CHR="$OPTARG";;
        h)      echo $USAGE
                exit 1;;
        esac
done

if [[ $VCF_FILE == "" || $EXON_FILE == "" || $CHR == "" ]]
then
        echo $USAGE
        exit 1
fi

DATAPATH=`dirname "$VCF_FILE"`
grep $CHECK1 $VCF_FILE | grep $CHECK2 | awk -v var1=$CHR '$1 == var1 {print $1 " " $2}' > $DATAPATH/hap.list.$CHR
cat $EXON_FILE | awk -v var1=$CHR '$1 == var1 {print $0}' > $DATAPATH/exon.list.$CHR

if [[ $? != 0 || ! -s $DATAPATH/hap.list.$CHR || ! -s $DATAPATH/exon.list.$CHR ]]
then
	echo SNV Count failed
	exit 1
fi

FD1=7
FD2=8
eof1=0
eof2=0

exec 7<$DATAPATH/hap.list.$CHR
exec 8<$DATAPATH/exon.list.$CHR
rm -f $DATAPATH/SNV.list.$CHR

SNV=-1
LOW=0
HIGH=0

while [[ $eof1 -eq 0 && $eof2 -eq 0 ]]
do
	while [[ $SNV -lt $LOW ]]
	do
		if read data1 <&$FD1
		then
			SNV=`echo "$data1" | cut -d ' ' -f 2`
		else
			eof1=1
			break
		fi
	done

	if [[ $SNV -le $HIGH ]]
	then
		echo SNV = $SNV TARGET REGION = $data2 >> $DATAPATH/SNV.list.$CHR
		if read data1 <&$FD1
		then
			SNV=`echo "$data1" | cut -d ' ' -f 2`
		else
			eof1=1
		fi
	else
		while [[ $SNV -gt $HIGH ]]
		do
			if read data2 <&$FD2
			then
                		LOW=`echo "$data2" | cut -f 2`
       	        		HIGH=`echo "$data2" | cut -f 3`
			else
				eof2=1
				break
			fi
		done
	fi
done

