#!/bin/bash

## global vars
proj=b2011142
mail="judith.felten@slu.se"

## stop on error
set -e

## check the arguments                                                                            
if [ $# != 1 ]; then
   echo "This script takes one argument: the tool which output is to be compressed"
   exit 1
fi


in=
case "$1" in
    STAR) 
	in=/proj/$proj/nobackup/STAR
	pattern="*Unmapped.out.mate[1-2]";;
esac

## test
if [ -z $in ]; then
    echo "The argument should be one of STAR or ...  [to be extended]"
    exit 1
fi  

## run
for f in `find $in -name $pattern -type f`
do 
line=`basename ${f//Unmapped.out.mate/}`
sbatch -A $proj --mail-user $mail -e $f.err -o $f.out -J Gzip-$line ../../../pipeline/runGzip.sh $f
done

