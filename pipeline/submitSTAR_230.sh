#!/bin/bash

## global vars
proj=b2011142
mail="judith.felten@slu.se"
intron=11000
gff=/proj/$proj/nobackup/reference/Populus-trichocarpa/v3.0/gff3/Ptrichocarpa_210_gene_exons.gff3
genome=/proj/$proj/nobackup/reference/Populus-trichocarpa/v3.0/indices/STAR

## two of the files need really long
time="06:00:00"


## stop on error
set -e

in=
case "$1" in
    sortmerna) 
	in=/proj/$proj/nobackup/sortmerna/230
	pattern="*_sortmerna_[1,2].fq.gz";;
    trimmomatic)
	in=/proj/$proj/nobackup/trimmomatic/230
	pattern="*_trimmomatic_[1,2].fq.gz";;
esac

## test
if [ -z $in ]; then
    echo "The argument should be one of sortmerna or trimmomatic"
    exit 1
fi  

## create the in and out dir
in=/proj/$proj/nobackup/$1/230
out=/proj/$proj/nobackup/STAR/$1/230
mkdir -p $out



# run                                                                                                                             
for f in `find $in -name "*.fq.gz" -type f`; do echo `basename ${f//_[1,2].fq.gz/}` ; done | sort | uniq | while read line;
do
sbatch -A b2012243  --mail-user $mail -e $out/$line.err -o $out/$line.out -J STAR-$line ../../../pipeline/runSTAR.sh -o $out -m $intro\
n $in/${line}_1.fq.gz $in/${line}_2.fq.gz $genome $gff -- --outQSconversionAdd -31 --outReadsUnmapped Fastx
done


