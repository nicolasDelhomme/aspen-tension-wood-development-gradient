#!/bin/bash

## stop on error
set -e

## args
mail="judith.felten@slu.se"
proj=b2011142
  

## create the in and out dir
in=/proj/$proj/nobackup/sortmerna/230
out=/proj/$proj/nobackup/trimmomatic/230
mkdir -p $out

## select all files
## NOTE THAT QUAL IS +33!!!! 
for f in `find $in -name "*.f*q.gz"`; do echo `basename ${f//_[1,2].f*q.gz*/}` ; done | sort | uniq | while read line;
do sbatch -A $proj --mail-user $mail -e $out/$line.err -o $out/$line.out -J Trim-$line ../../../pipeline/runTrimmomatic.sh $in/${line}_1.f*q.gz $in/${line}_2.f*q.gz $out 33 SLIDINGWINDOW:5:20 MINLEN:50
done

