#!/bin/bash
#$ -cwd


while getopts l:s:g:p:h opt
  do 
  case "$opt" in
      l) list="$OPTARG";;
      s) snp="$OPTARG";;
      g) geno="$OPTARG";;
      p) pheno="$OPTARG";;
      h) echo $USAGE
        exit 1;;
  esac
done

OUTDIR=$list"_VT-test"
mkdir -p $OUTDIR

for f in `cat $list`
  do 
  egrep -w "^$f" $snp  | cut -f2,3 > $OUTDIR/vcf.$f.csnp; 
  egrep -w "^$f" $geno | cut -f2,3,4 > $OUTDIR/vcf.$f.cgeno
  
  
  Rscript --slave /ifs/home/c2b2/af_lab/saec/code/ExomePipe/association/rareVariantTests.R -p 1000000 -a $pheno -b $OUTDIR/vcf.$f.csnp -c $OUTDIR/vcf.$f.cgeno  --multicore  > $OUTDIR/vt.$f.R-result
done
