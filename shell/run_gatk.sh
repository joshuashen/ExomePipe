PRE='pilot'
SRC='src'

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

# Merge SNP calls
OUT="$PRE.calls"
qsub gatk_varcalibrate.scr $OUT/$PRE
