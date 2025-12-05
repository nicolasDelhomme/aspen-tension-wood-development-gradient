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
	in=/proj/$proj/nobackup/sortmerna
	pattern="*_sortmerna_[1,2].fq.gz";;
    trimmomatic)
	in=/proj/$proj/nobackup/trimmomatic
	pattern="*_trimmomatic_[1,2].fq.gz";;
esac

## test
if [ -z $in ]; then
    echo "The argument should be one of sortmerna or trimmomatic"
    exit 1
fi  

## create the in and out dir
out=/proj/$proj/nobackup/STAR/$1
mkdir -p $out

## run
for d in `find $in -maxdepth 1 -mindepth 1 -type d`
do
mkdir -p $out/`basename $d`
## then per file
for f in `find $d -name $pattern -type f`; do echo `basename ${f//_[1,2].fq.gz/}` ; done | sort | uniq | while read line;
do 
sbatch -A $proj -t $time --mail-user $mail -e $out/$line.err -o $out/$line.out -J STAR-$line ../../../pipeline/runSTAR.sh -o $out/`basename $d` -m $intron $d/${line}_1.fq.gz $d/${line}_2.fq.gz $genome $gff -- --outQSconversionAdd -31 --outReadsUnmapped Fastx
done
done

