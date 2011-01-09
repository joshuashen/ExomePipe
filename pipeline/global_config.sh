# Global config file for some default parameters

# JAVA related parameters
export HEAP=4000
export JAVA="java -Xmx${HEAP}m -Djava.io.tmpdir=/tmp"

# Locations for all tools
export GATK="$JAVA -jar /ifs/data/c2b2/ip_lab/shares/SOFTWARE/Sting/dist/GenomeAnalysisTK.jar"
export FIXMATE="$JAVA -jar /ifs/data/c2b2/ip_lab/shares/SOFTWARE/picard-tools-1.29/FixMateInformation.jar"
export SAMTOOLS="/ifs/home/c2b2/ip_lab/yshen/usr/bin/samtools"

# This contains the base path to where all the scripts are available. Modify this as per the clone location..
export BPATH="/ifs/scratch/c2b2/ip_lab/aps2157/ExomePipe/pipeline/"

# Do not modify the variables below.
export ERRORMESSAGE="#### ERROR"
export ERRORMESSAGE1="The following error has occurred"
export EXONUSAGE="Please specify either the file containing the interval list using -E or the sequences using -L"
export JOB_STAT="Following"

# Placeholders. Modify the ExonList variable if needed.
export INP=""
export CHR=""
export REF=""
export DBSNP=""
export Platform=""
export List=""
export ExonList="/ifs/scratch/c2b2/ip_lab/aps2157/ExomePipe/pipeline/data/old_samples/captureexoncoordinatesonly.txt"

# Do not modify the variables below unless locations change.
export REF_V37="/ifs/data/c2b2/ip_lab/shares/DATA/Sequencing/resources/wu_build36.fasta"
export DBSNP_V37="/ifs/scratch/c2b2/ip_lab/aps2157/ExomePipe/pipeline/data/TCGA_WU/dbsnp_129_b36.ordered.rod"
export REF_HG18="/ifs/data/c2b2/ip_lab/shares/DATA/Sequencing/resources/bcm_hg18.fasta"
export DBSNP_HG18="/ifs/scratch/c2b2/ip_lab/aps2157/ExomePipe/pipeline/data/old_samples/dbsnp_129_hg18.ordered.rod"

# Only make additions below.

export TCGA_REF="${REF_V37}"
export TCGA_DBSNP="${DBSNP_V37}"
export TCGA_Platform="illumina"
export TCGA_List="${ExonList}"

export YALE_REF="${REF_HG18}"
export YALE_DBSNP="${DBSNP_HG18}"
export YALE_Platform="illumina"
export YALE_List="${ExonList}"

