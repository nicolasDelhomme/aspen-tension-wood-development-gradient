#!/bin/bash -l

set -ex

proj=u2015018
mail=carolin.seyfferth@slu.se
#mail=nicolas.delhomme@umu.se
in=/mnt/picea/projects/aspseq/jfelten/06_ERF_Project/raw
out=/mnt/picea/projects/aspseq/jfelten/06_ERF_Project
genome=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/STAR/2.5.2b/Potra01
gff3=/mnt/picea/storage/reference/Populus-tremula/v1.1/gff3/Potra01-gene-synthetic-transcripts-wo-intron.gff3
gtf=/mnt/picea/storage/reference/Populus-tremula/v1.1/gtf/Potra01-gene-mRNA-wo-intron.gtf
kallisto_fasta=/mnt/picea/storage/reference/Populus-tremula/v1.1/fasta/Potra01-mRNA.fa
kallisto_index=/mnt/picea/storage/reference/Populus-tremula/v1.1/indices/kallisto/Potra01-mRNA.fa.inx
start=4
end=9
mem=100G

module load bioinfo-tools FastQC sortmerna Trimmomatic star htseq samtools kallisto

if [ -z $UPSCb ]; then
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir."
fi

for line in ( 3_111130_BD07JJACXX_t8-4_index4 4_110707_BB0B0HABXX_k5-1_ind_4 5_110707_BB0B0HABXX_k8-3_ind_2 5_110707_BB0B0HABXX_t8-2_ind_5 5_110707_BB0B0HABXX_t8-3_ind_6 6_110707_BB0B0HABXX_t8-4_ind_4 6_110707_BB0B0HABXX_t8-5_ind_5 7_110707_BB0B0HABXX_t11-4_ind_4 7_110707_BB0B0HABXX_t11-4_ind_4 8_110815_AB0B10ABXX_t11-4_ind_4 ); do 
  bash $UPSCb/pipeline/runRNASeqPreprocessing.sh -s $start -e $end -g $genome \
  -f $kallisto_fasta -K $kallisto_index \
  -G $gtf -H $gff3 -m $mem $proj $mail $in/${line}_1.fastq.gz $in/${line}_2.fastq.gz $out
done
