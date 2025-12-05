#!/bin/bash

## global vars
proj=b2011142
mail="judith.felten@slu.se"
gff=/proj/$proj/nobackup/reference/Populus-trichocarpa/v3.0/gff3/Ptrichocarpa_210_gene_exons.min_1000.const_exons.gff
in=/proj/$proj/nobackup/STAR/sortmerna

## stop on error
set -e

## create the in and out dir
out=/proj/$proj/nobackup/MISO/PE_distribution
mkdir -p $out

## run
for d in `find $in -maxdepth 1 -mindepth 1 -type d`
do
mkdir -p $out/`basename $d`
## then per file
for f in `find $d -name "*.bam" -type f`
do 
line=`basename ${f//.bam/}`
sbatch -A $proj --mail-user $mail -e $out/$line.err -o $out/$line.out -J MISO-PE-$line ../../../pipeline/runMisoPeDistribution.sh $out/`basename $d` $f $gff
done
done

