PRE='pilot'
SRC='src'

SAMPLES='dot_20090601_1 dot_20090601_2 dot_20091217_1 dot_20091217_2 ellen_20100209_1 ellen_20100209_2 helen_20090605_1 helen_20090605_2 helen_20100108_1 helen_20100108_2 tomoe_20090714_1'

# Base quality score recalibration
OUT="$PRE.recalibrated"
mkdir $OUT
for s in $SAMPLES; do
  ls --color=never $SRC/$s.*bam > $OUT/$s.list
  qsub gatk_calibrate.scr $OUT/$s
done
exit

OUT="$PRE.recalibrated"
for s in $SAMPLES; do
  ls --color=never $SRC/$s.*bam > $OUT/$s.list
  for c in `seq 1 22`; do
    qsub gatk_recalibrate.scr $OUT/$s $c
  done
done
exit

# Local realignment around indels
OUT="$PRE.realigned"
PREV="$PRE.recalibrated"
mkdir $OUT
for c in `seq 1 22`; do
  ls --color=never $PREV/*.$c.recalibrated.bam > $OUT/$PRE.$c.list
  qsub gatk_realign.scr $OUT/$PRE.$c $c
done
exit

# Indel & SNP calling
OUT="$PRE.calls"
PREV="$PRE.realigned"
mkdir $OUT
for c in `seq 1 22`; do
  echo $PREV/$PRE.$c.fixed.bam > $OUT/$PRE.$c.list
  # Indel calls
  qsub gatk_indelcall.scr $OUT/$PRE.$c $c
  # SNP calls
  qsub gatk_snpcall.scr $OUT/$PRE.$c $c
done
exit

# Filter SNP calls
OUT="$PRE.calls"
for c in `seq 1 22`; do
  qsub gatk_filter.scr $OUT/$PRE.$c $c
done
exit

# Merge and recalibrate variant quality score
OUT="$PRE.calls"
qsub gatk_merge.scr "$OUT/$PRE.*.snps.filtered.vcf" $OUT/$PRE.snps.filtered
exit
qsub gatk_varcalibrate.scr $OUT/$PRE.snps.filtered
exit

# Impute variants with BEAGLE
PREV="$PRE.calls"
OUT="$PRE.beagle"
mkdir $OUT
ln -s $PREV/$PRE.snps.filtered.recalibrator_output.filtered.vcf* $OUT
for c in `seq 1 22`; do
  qsub gatk_beagle.scr $OUT/$PRE.recalibrator.output.filtered $c
done
