#!/bin/bash                                                                                                                           
## vars                                                                                                                               
proj=b2011142
mail="judith.felten@slu.se"

## check the arguments                                                                                                                
if [ $# != 1 ]; then
   echo "This script takes one argument: the input data  name"
   exit 1
fi

in=
nam="*.fq.gz"
typ="f"
case "$1" in
    raw)
	in=/proj/$proj/sequence_data/raw/soft_linked/230
        typ="l"
        nam="*.fastq.gz";;
    trimmomatic)
	in=/proj/$proj/nobackup/trimmomatic/230
        typ="f"
        nam="*.fq.gz";;
    sortmerna)
	in=/proj/$proj/nobackup/sortmerna/230
        typ="f"
        nam="*.fq.gz";;
esac

if [ -z $in ]; then
    echo "The argument should be one of raw, trimmomatic or sortmerna"
    exit 1
fi

## check $in                                                                                                                          
if [ ! -d $in ]; then
    echo "The tool directory $1 does not exist. Make sure that the tool was run first."
    exit 1
fi

## create the dir                                                                                                                     
out=/proj/$proj/qa/$1select 
mkdir -p $out

for file in `find $in -name $nam -type $typ`; do
fnam=`basename $file`
fnam=${fnam//.f*q.gz/}
sbatch -A $proj --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J FastQC-$fnam ../../../pipeline/runFastQC.sh $out $file
done