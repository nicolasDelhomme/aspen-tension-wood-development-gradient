## Judith's transcript

##copy files into R (aspseq) - only for the first time the analysis is done
##make folders under aspseq and go there, then
## scp -r judith@kalkyl.uppmax.uu.se:/proj/b2012243/nobackup/data/tensionwood/DESeq_Count/*.log .


##then in bash do on the folder on aspseq
## cd ../accdata_SeqCount/
## judith@aspseq:~/ethylene-insensitive/accdata_SeqCount$ cat ./*.log | awk '{ print $1 }' | sort -u > accdata_gene_list.log


## R
## Create the big table with all expression data
## acctable<-read.table(file="ethylene-insensitive//accdata_SeqCount//accdata_gene_list.log")
## head(acctable)

## ACC_T89_3_table<-read.table(file="ethylene-insensitive//accdata_SeqCount//T89-10h-ACC3_trimmomatic_sortmerna_STAR_STAR.bam_counts.log")
## for all files

## Then create the big table by associating the respective count files to the acctable
## acctable$ACC_1E_2<- ACC_1E_2_table$V2[match(acctable$V1,ACC_1E_2_table$V1)]
## acctable$ACC_1E_3<- ACC_1E_3_table$V2[match(acctable$V1,ACC_1E_3_table$V1)]

## count_table <- data.frame(acctable[1:66223,2:19])
## row.names(count_table) <- acctable$V1
## head(count_table)

##Define sample classes/conditions
## meta<-data.frame(row.names=colnames(count_table),genotype=c(rep("1E",3),rep("T89",3),rep("X6",3),rep("1E",3),rep("T89",3),rep("X6",3)),condition=c(rep("NW",9),rep("TW",9)))
## meta<-data.frame(row.names=colnames(count_table),condition=c(rep("ACC",3),rep("H2O",3),rep("ACC",3),rep("H2O",3),rep("ACC",3),rep("H2O",3)),genotype=c(rep("1E",6),rep("X6",6),rep("T89",6)))


##Or with Nico's script:
## library
library(RColorBrewer)

## set the work dir to wherever the data is
setwd("~/ethylene-insensitive-HTSeq_DESeq/accdata/inputdata")

## read the HTSeq files and create the count table
res <- lapply(dir(".",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
})

## name the samples
names(res) <- sub("_trimmomatic_sortmerna_STAR.txt","",dir(".",pattern="*.txt"))
## SANITY check to ensure the gene order is similar in all files
stopifnot(all(sapply(res[2:length(res)],function(r,o){
  all(r[,1] == o)
},res[[1]][,1])))
## extract the count.table
addInfo <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- sub("\\.0$","",res[[1]][,1])[-sel]
##count.table <- count.table[grep("Potri",rownames(count.table)),]

### ==============================
## get the last stat lines
### ==============================
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- addInfo
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]
count.stats
pal=brewer.pal(6,"Dark2")[1:nrow(count.stats)]
barplot(as.matrix(count.stats),beside=TRUE,las=2,
        main="read proportion",
        ylim=c(1,2e7),col=pal)
legend("top",fill=pal,legend=rownames(count.stats),bty="n",cex=0.5,horiz=TRUE)
##Then we continue like before NOTE: In Nico's script the count.table is what with Bastian we named count_table
## Define sample classes/conditions
## for TW data
meta<-data.frame(row.names=colnames(count_table),genotype=c(rep("1E",3),rep("T89",3),rep("X6",3),rep("1E",3),rep("T89",3),rep("X6",3)),condition=c(rep("NW",9),rep("TW",9)))

##for ACC data
meta<-data.frame(row.names=colnames(count.table),condition=c(rep("ACC",3),rep("H2O",3),rep("ACC",3),rep("H2O",3),rep("ACC",3),rep("H2O",3)),genotype=c(rep("1E",6),rep("T89",6),rep("X6",6)))

## to check that its 1EACC, 1EACC, 1EACC, 1EH20,… and that there's no effect on normalizing only two conditions together
#doesnt work for whatever reason conditions <- paste(c(rep("1E",6),rep("X6",6),rep("T89”,6)),c(rep("ACC",3),rep("H2O",3),rep("ACC",3),rep("H2O",3),rep("ACC",3),rep("H2O”,3)),sep=“”)
conditions<-data.frame(row.names=colnames(count.table),condition=c(rep("1EACC",3),rep("1EH2O",3),rep("T89ACC",3),rep("T89H2O",3),rep("X6ACC",3),rep("X6H2O",3)))

library("DESeq")

##Quality check on all data incl. normalization
#count_table[is.na(count_table)] <- 0 
## count.table[is.na(count.table)]<- 0
sum(rowSums(is.na(count.table)) > 0)
cdsFull<-newCountDataSet(count.table,as.character(conditions$condition))
cdsFull <- estimateSizeFactors(cdsFull)
sizeFactors(cdsFull)
cdsFullBlind <- estimateDispersions(cdsFull, method="blind")
cdsFull <- estimateDispersions(cdsFull)

vsdFull <- varianceStabilizingTransformation(cdsFullBlind)
library("RColorBrewer")
library("gplots")

#chose where results shall be saved from now on
setwd("~/ethylene-insensitive-HTSeq_DESeq/accdata/Overall data comparison_Feb2014/")

select = order(rowMeans(counts(cdsFull)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
png(filename="Heatmap_top30_vst_accdata.png")
heatmap.2(exprs(vsdFull)[select,], col = hmcol, trace="none", margin=c(10, 6))
dev.off()

##Sample to Sample distance
dists = dist( t( exprs(vsdFull) ) )
mat = as.matrix( dists )
rownames(mat) = colnames(mat) = with(pData(cdsFullBlind), paste(conditions))

postscript(file="sample-to-sample-distances_accdata.eps", width=25.6,height=14.4)
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
dev.off()

##PCA plot
cdsM<-newCountDataSet(count.table,meta)
cdsM <- estimateSizeFactors(cdsM)
sizeFactors(cdsM)
cdsM <- estimateDispersions(cdsM, method="blind")
vsdM <- varianceStabilizingTransformation(cdsM)
postscript(file="PCA_plot_ACC-vs-H2O.eps", width=25.6,height=16.4)
print(plotPCA(vsdM, intgroup=c("condition", "genotype")))
dev.off()

##go to sample wise analysis 
##For dispersion and differential gene expression analysis for treatment differences
## samples<-conditions$conditions
## countTable<-count.table[,samples]
#condition=conditions$conditions
## cds<-newCountDataSet(count.table, conditions$conditions)
## cds<- estimateSizeFactors(cdsFull)
## sizeFactors(cds)
## cds = estimateDispersions(cds)
png(filename="dispersion_estimation_1E-ACC-H2O.png", width=1920,height=1080,units="px",pointsize=20)
plotDispEsts(cdsFull)
dev.off()

source("~/UPSCb/src/R/plotDispLSD.R")
plotDispLSD(cdsFull)

##Differential expression and Histogram on log2 and p-Values
res = nbinomTest( cdsFull, "1EH2O", "1EACC" )
png(filename="MAplot_1E-ACC-vs-H2O.png", width=1920,height=1080,units="px",pointsize=20)
plotMA(res)
dev.off()

postscript(file="histogram-pvalues_1E-ACC-vs-H2O.eps", width=25.6,height=14.4)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

##Save the data for each comparison (genotype)
write.csv(res,file="res_1E-ACC-vs-H2O.csv")

res = nbinomTest( cdsFull, "T89H2O", "T89ACC" )
plot(res$baseMean,res$log2FoldChange,log="x",pch=ifelse(res$padj<.05,19,20),col=ifelse(res$padj<.05,"red","black"))
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
plotMA(res)

res = nbinomTest( cdsFull, "T89ACC", "1EACC" )
plotMA(res)



head(count.table)
cdsT <- newCountDataSet(count.table[,c(1:6,13:18)],conditions$condition[c(1:6,13:18)])
cdsT <- estimateSizeFactors(cdsT)
cdsT <- estimateDispersions(cdsT)
cdsT
sampleNames(cdsT)
conditions(cdsT)
plotDispLSD(cdsT)
res = nbinomTest( cdsT, "1EH2O", "1EACC" )
plotMA(res)



library(DESeq)
cds1E <- newCountDataSet(count.table[,c(1:6)],conditions$condition[c(1:6)])
cds1E <- estimateSizeFactors(cds1E)
cds1E <- estimateDispersions(cds1E)
res = nbinomTest( cds1E, "1EH2O", "1EACC" )
plotMA(res)

str(res)

plot(res$baseMean,res$log2FoldChange,log="x",pch=ifelse(res$padj<.05,19,20),col=ifelse(res$padj<.05,"red","black"))



---
  ######For genotype comparison
  
  ##For dispersion and differential gene expression analysis for genotype differences
  samples<-meta$condition =="H2O"
countTable<-count.table[,samples]
genotypes=meta$genotype[samples]
cds<-newCountDataSet(countTable, genotypes)
cds<- estimateSizeFactors(cds)
sizeFactors(cds)
cds = estimateDispersions( cds )
png(filename="dispersion_estimation_H2O-T89-X6-1E.png", width=1920,height=1080,units="px",pointsize=20)
plotDispEsts(cds)
dev.off()

##Differential expression and Histogram on log2 and p-Values
res = nbinomTest( cds, "X6", "1E" )
png(filename="MAplot_H2O-1E-vs-H2O-X6.png", width=1920,height=1080,units="px",pointsize=20)
plotMA(res)  
dev.off()

postscript(file="histogram-pvalues_H2O-1E-vs-H2O-X6.eps", width=25.6,height=14.4)
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
dev.off()

##Save the data for each comparison (genotype)
write.csv(res,file="res_H2O-1E-vs-H2O-X6.csv")

----
  
  ##Save the count table (only once because all data is in there)
  write.csv(count.table,file="count.table_TW_data.csv")	

##in the end copy all eps and csv and png files to your computer
scp judith@aspseq.fysbot.umu.se:ethylene-insensitive/tensionwood_SeqCount/*.eps .


##To match CSV list to annotation file

geneannotation<-read.delim(file="ethylene-insensitive/Ptrichocarpa_210_annotation_info.txt", header=FALSE)

(HArdy's script, didn't work for me)
T89results.annots<- NULL
for (i in seq(1,length(T89results[,1]))){
  num=match(T89results$id[i],geneannotation$V4)
  geneannotation.tmp=geneannotation[num,4:13]
  T89results.annots=rbind(T89results.annots,geneannotation.tmp)
}
T89results.annots=cbind(T89results.annots,T89results) 


geneannotation$id_Results<- count.table$rownames[match(geneannotation$V2,count.table$rownames)]