#!/bin/bash                                                                                                                                                                                                              

proj=b2011142
tmp=/glob/$USER/tmp
mkdir -p $tmp
data=/glob/$USER/data
mail="judith.felten@slu.se"
## one of the file takes really long                                                                                                                                                                                          
time="24:00:00"

## error                                                                                                                                                                                                                 
set -e


## this script takes on argument                                                                                                                                                                                              
in=
case "$1" in
    raw)
        in=/proj/$proj/sequence_data/raw/qa_passed
        typ=l;;
    trimmomatic) in= /proj/$proj/sequence_data/nobackup/trimmomatic
        typ=f;;
esac

## test                                                                                                                                                                                                                       
if [ -z $in ]; then
    echo "The argument should be one of raw or trimmomatic"
    exit 1
fi

## create the in and out dir                                                                                                                                                                                                  
out=/proj/$proj/nobackup/sortmerna
mkdir -p $out

## check that the needed data is in the user glob                                                                                                                                                                             
if [ ! -e ../../../data/sortmerna ]; then
    if [ ! -d $data/sortmerna ]; then
        mkdir -p $data/sortmerna
        echo -n "Enter your username at picea.plantphys.umu.se: "
        read picea
        scp -P 922 -r $picea@picea.plantphys.umu.se:/mnt/picea/storage/reference/rRNA/sortmerna $data/sortmerna
    fi
    ln -sf $data/sortmerna -T ../../../data/sortmerna
else
    if [ ! -L ../../../data/sortmerna ]; then
        echo "Aborting as we do not want to overwrite data." `pwd`"/../../../data/sortmerna should be a symbolic link."
        exit 1
    fi
fi

## then start the script                                                                                                                                                                                                      
## loop per dir                                                                                                                                                                                                               
for d in `find $in -maxdepth 1 -mindepth 1 -type d`
do
mkdir -p $out/`basename $d`
## then per file                                                                                                                                                                                                              
for f in `find $d -name "*_[1,2].fastq.gz" -type $typ`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
fnam=`basename $line`
sbatch -A $proj -t $time --mail-user $mail -e $out/$fnam.err -o $out/$fnam.out -J smr-$fnam ../../../pipeline/runSortmerna.sh $out/`basename $d` $tmp ${line}_1.fastq.gz ${line}_2.fastq.gz
done
done


