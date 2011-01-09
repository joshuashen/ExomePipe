###set global config
export REFTYPE="build"  ## build: no "chr" in chromosome names, aka >1 >2 etc;  hg: >chr1 >chr2 etc.
export REF="/ifs/data/c2b2/af_lab/saec/Sequencing/resources/references/wu_build36.fasta"
export ExonFile="/ifs/data/c2b2/af_lab/saec/Sequencing/resources/exomes/captureexoncoordinatesonly.txt"
export DBSNP="/ifs/data/c2b2/af_lab/saec/Sequencing/resources/references/dbsnp_130_b36.rod.wu" 
export SAMTOOLS="/ifs/home/c2b2/ip_lab/yshen/usr/bin/samtools"
###  export TEMP="/ifs/scratch/c2b2/af_lab/saec/temp/"  # temp dir for Java
export GATKJAR="/ifs/data/c2b2/af_lab/saec/Software/Sting/dist/GenomeAnalysisTK.jar"
export PICARD="/ifs/data/c2b2/af_lab/saec/Software/picard-tools-1.33"
export BWA="/ifs/home/c2b2/ip_lab/yshen/usr/bin/bwa"  # bwa 0.5.8
export BPATH="/ifs/home/c2b2/af_lab/saec/code/ExomePipe/calling/"
export STING="/ifs/data/c2b2/af_lab/saec/Software/Sting"
export HEAP=4000  # default 4G heap 

