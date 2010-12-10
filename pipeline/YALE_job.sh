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

while getopts I:h o
do      case "$o" in
        I)      INP="$OPTARG";;
        h)      echo $USAGE
                exit 1;;
        esac
done

if [[ $INP == "" ]]
then
	echo "Please pass input bam file as parameter"
	exit 1
fi

DATAPATH=`dirname "$INP"`
rm -f $DATAPATH/*.output

$BPATH/run_pipeline.sh -I $INP -R $YALE_REF -E $YALE_List -D $YALE_DBSNP -P $YALE_Platform > $DATAPATH/pipeline.output

