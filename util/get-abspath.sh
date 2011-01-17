#!/bin/bash
#$ -cwd

in=$1

for f in `cat $in`
  do
  g=`readlink -e $f`
  echo $g
done
