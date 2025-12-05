#!/bin/bash -l

set -ex

proj=u2015018
#mail=carolin.seyfferth@slu.se
mail=nicolas.delhomme@umu.se
in=/mnt/picea/projects/aspseq/jfelten/06_ERF_Project/raw
out=/mnt/picea/projects/aspseq/jfelten/06_ERF_Project
genome=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/STAR/2.5.2b/Potra01
gff3=/mnt/picea/storage/reference/Populus-tremula/v1.1/gff3/Potra01-gene-synthetic-transcripts-wo-intron.gff3
gtf=/mnt/picea/storage/reference/Populus-tremula/v1.1/gtf/Potra01-gene-mRNA-wo-intron.gtf
start=7
end=9
mem=100G

module load bioinfo-tools FastQC sortmerna Trimmomatic star htseq samtools kallisto

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for f in `find $in -name "*_[1,2].fastq.gz"`; do echo "${f//_[1,2].fastq.gz/}" ; done | sort | uniq | while read line;
do
  bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -g $genome \
  -G $gtf -H $gff3 -m $mem $proj $mail ${line}_1.fastq.gz ${line}_2.fastq.gz $out
done
