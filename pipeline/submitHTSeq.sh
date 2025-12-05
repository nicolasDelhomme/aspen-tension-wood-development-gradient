#!/bin/bash

proj=b2011142
tmp=/glob/$USER/tmp
gff3=/proj/$proj/nobackup/reference/Populus-trichocarpa/v3.0/gff3/Ptrichocarpa_v3.0_210_synthetic-gene-models.gff3
mail="judith.felten@slu.se"

## some take long
tim="24:00:00"

## stop on error
set -e

## what index
in=
out=
case "$1" in
    erf) 
	in=/proj/$proj/data/ERF/STAR/trimmomatic/soft-linked
	out=/proj/$proj/analysis/ERF/HTSeq
	job="ERF";;
    tensionwood) 
	in=/proj/$proj/data/TW_sections_MK/STAR/trimmomatic/soft-linked
	out=/proj/$proj/analysis/TW_sections_MK/HTSeq
	job="TW-MK";;
esac

if [ -z $in ]; then
    echo "The argument should be one of erf or tensionwood"
    exit 1
fi

if [ ! -d $out ]; then
    mkdir -p $out
fi

if [ ! -d $tmp ]; then
    mkdir -p $tmp
fi

for f in `find $in -name "*.bam" -type l`
do 
fnam=`basename ${f//.bam/}`
sbatch -A $proj -t $tim --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J htseq-$job-$fnam ../../../pipeline/runHTSeq.sh $out $tmp $f $gff3
done

