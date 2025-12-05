#!/bin/bash

## global vars
proj=b2011142
mail="judith.felten@slu.se"
inx=/proj/$proj/nobackup/reference/Populus-trichocarpa/v3.0/indices/MISO
in=/proj/$proj/nobackup/STAR/sortmerna

## stop on error
set -e

## create the in and out dir
out=/proj/$proj/nobackup/MISO
mkdir -p $out

## run
for d in `find $in -maxdepth 1 -mindepth 1 -type d`
do
    mkdir -p $out/`basename $d`
## then per file
    for f in `find $d -name "*M_*.bam" -type f`
    do 
	line=`basename ${f//.bam/}`
	mkdir -p $out/`basename $d`/$line
	sbatch -A $proj --mail-user $mail -e $out/$line.err -o $out/$line.out -J MISO-$line ../../../pipeline/runMiso.sh $out/`basename $d`/$line $inx $f -l 101 -m 175 -s 50
    done
done

