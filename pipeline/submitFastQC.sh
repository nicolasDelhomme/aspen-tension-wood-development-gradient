#!/bin/bash

echo "Change me to use the Git runner!\n"
exit 1;

## check $1
od=
case "$1" in
    raw) od="/proj/b2011142/sequence_data/raw/soft_linked";;
    trimmed) od="/proj/b2011142/sequence_data/nobackup/trimmed";;
esac

if [ -z $od ]; then
    echo "The first argument should be one of 'raw' or  'trimmed'"
    exit 1
fi

for d in `find $od -maxdepth 1 -mindepth 1 -type d`
# for d in `dir`
do
if [ "$1" == "raw" ]; then
    for f in `ls $d`; do echo "${f//_[1,2].fastq.gz/}" ; done | uniq | while read line;
    do sbatch -e /proj/b2011142/qa/$1/$line.err -o /proj/b2011142/qa/$1/$line.out -J QA-$line /proj/b2011142/pipeline/runFastQC.sh $d/$line* $1
    done
else
    for f in `ls $d`; do echo "${f//_[1,2]_*.fq.gz/}" ; done | uniq | while read line;
    do sbatch -e /proj/b2011142/qa/$1/$line.err -o /proj/b2011142/qa/$1/$line.out -J QA-$line /proj/b2011142/pipeline/runFastQC.sh $d/$line* $1
    done
fi
done
