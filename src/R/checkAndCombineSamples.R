## in R
library(HistoneChIPseq)
library(LSD)

## setwd
setwd("/mnt/picea/storage/projects/06_ERF_Project")

## read the bam files
bamFiles <-BamFileList(dir("data/STAR/sortmerna/",pattern="*.bam$",full.names=TRUE),dir("data/STAR/sortmerna/",pattern="*.bam$",full.names=TRUE))

## sum in 25bp bins
sexp <- HistoneChIPseq(bamFiles=bamFiles,
                       binSize=25L,
                       chr.sel=sprintf("Chr%02d",1:19),
                       method="shift",
                       shift=0L,
                       strand="both",
                       yieldSize=10^6L,
                       nnodes=12L,
                       normalize=TRUE)

## perform the comparison: 14 comparisons
smpl <- factor(gsub("[2-8]_[0-9]+_[A-Z0-9]{10}_|_ind.*","",colData(sexp)$FileName))

## work per sample
sapply(levels(smpl),function(smp,smpl,sexp){
  apply(combn(which(smpl==smp),2),2,function(co,sexp){

    x <- log10(unlist(assays(sexp)$RPGC[,co[1]][[1]]))
    y <- log10(unlist(assays(sexp)$RPGC[,co[2]][[1]]))
    xy.min <- min(c(x[!is.infinite(x)],y[!is.infinite(y)]))
    x[is.infinite(x)] <- xy.min -1
    y[is.infinite(y)] <- xy.min -1
    png(file.path("analysis","replicate-comparison",paste(colData(sexp)$FileName[co[1]],"-",colData(sexp)$FileName[co[2]],".png",sep="")),width=600,height=600,pointsize=16)
    comparisonplot(x,y,ab=TRUE,
                   xlab=colData(sexp)$FileName[co[1]],
                   ylab=colData(sexp)$FileName[co[2]],
                   main=smp)
    dev.off()
  },sexp)},smpl,sexp)


