#!/bin/bash
#$ -cwd

perm1=10000 
perm2=100000000
pcut=0.05

dirbase="/ifs/home/c2b2/af_lab/saec/code/ExomePipe/association"

while getopts l:s:g:p:b:h opt
  do 
  case "$opt" in
      l) list="$OPTARG";;
      s) snp="$OPTARG";;
      g) geno="$OPTARG";;
      p) pheno="$OPTARG";;
      b) dirbase="$OPTARG";;
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
  
  # first do permuation $perm1 times
  ruby $dirbase/vt-test-engine.rb  -a $pheno -b $OUTDIR/vcf.$f.csnp -c $OUTDIR/vcf.$f.cgeno -o $OUTDIR/vt.$f.R-result -d $dirbase -g $f

#   Rscript --slave $dirbase/rareVariantTests.R -p $perm1  -a $pheno -b $OUTDIR/vcf.$f.csnp -c $OUTDIR/vcf.$f.cgeno  --multicore  > $OUTDIR/vt.$f.R-result

  
done
