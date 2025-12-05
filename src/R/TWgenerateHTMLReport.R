### ==============================
## load the necessary libraries
### ==============================
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(scatterplot3d))
source("~/UPSCb/src/R/plot.multidensity.R")
source("~/UPSCb/src/R/rmd.R")

### ==============================
## set the working dir
### ==============================
setwd("~/06_ERF_Project/")

### ==============================
## read the HTSeq files in a matrix
### ==============================
res <- mclapply(dir("analysis/HTSeq/tensionwood",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=9)
names(res) <- sub("\\.txt","",dir("analysis/HTSeq/tensionwood",pattern="*.txt"))

### ==============================
## get the count table
### ==============================
addInfo <- c("no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",2))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]

### ==============================
## get the last stat lines
### ==============================
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- addInfo
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

### ==============================
## plot the stats
### ==============================
pal=brewer.pal(6,"Dark2")[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=pal,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,2e+6))
legend("top",fill=pal,legend=rownames(count.stats),bty="n",cex=0.8)

### ==============================
## 21% of the genes are not expressed
## from a total of 41,335 genes
### ==============================
sel <- rowSums(count.table) == 0
sum(sel) / nrow(count.table)
length(sel)

### ==============================
## display the per-gene mean expression
## i.e. the mean raw count of every 
## gene across samples is calculated
## and displayed on a log10 scale
### ==============================
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

### ==============================
## The same is done for the individual
## samples colored by sample type
### ==============================
pal=brewer.pal(4,"Dark2")
plot.multidensity(log10(count.table),col=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
                  legend.x="topright",legend.cex=0.5,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)",lwd=3)

### ==============================
## For visualization, the data is
## submitted to a variance stabilization
## transformation
### ==============================
conditions <- sub("-[1-6]","",colnames(count.table))
dds <- DESeqDataSetFromMatrix(countData = count.table,
                              colData = data.frame(condition=conditions),
                              design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

### ==============================
## what about the library sizes?
## they are all over the place...
## let's check the relative log transformation also
### ==============================
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

## compare rld and and vsd
rld <- rlogTransformation(dds, blind=TRUE)
px     <- counts(dds)[,1] / sizeFactors(dds)[1]
ord    <- order(px)
ord    <- ord[px[ord] < 150]
ord    <- ord[seq(1, length(ord), length=50)]
last   <- ord[length(ord)]
vstcol <- c("blue", "black")
matplot(px[ord],
        cbind(assay(vsd)[, 1], log2(px))[ord, ],
        type="l", lty=1, col=vstcol, xlab="n", ylab="f(n)")
legend("bottomright",
       legend = c(
         expression("variance stabilizing transformation"),
         expression(log[2](n/s[1]))),
       fill=vstcol)

### ==============================
## plot the mean agains the variance (sd)
## it should be flat if the VST worked properly
## it is not bu it is likely due to the 
## variability of the gene expression
## in normal wood compared to tension wood
## This is not a concern for visualization
## but might be for normalization
### ==============================
meanSdPlot(assay(vsd)[rowSums(count.table)>0,], ylim = c(0,2.5))

#plot.multidensity(as(assay(vsd),"list"))

###============================
## select a cutoff of 1 and plot the mean against the variance again.
## Cutoff 1 not  enough for getting an equal variance.
###============================

sel <- rowMeans(assay(vsd)) > 1
meanSdPlot(assay(vsd)[rowMeans(assay(vsd)) > 1,])

###============================
## select a cutoff of 2 and plot the mean against the variance again.
## looks pretty ok for getting and equal variance.
###============================

sel <- rowMeans(assay(vsd)) > 2
meanSdPlot(assay(vsd)[rowMeans(assay(vsd)) > 2,])

###============================
## select a cutoff of 5 and plot the mean against the variance again
## not so much improvement for getting equal variance.
## Cutoff 2 may be sufficient.
###============================

sel <- rowMeans(assay(vsd)) > 5
meanSdPlot(assay(vsd)[rowMeans(assay(vsd)) > 5,])

### ==============================
## Perform a Principal Component Analysis
### ==============================
pc <- prcomp(t(assay(vsd)))
percent <- round(summary(pc)$importance[2,]*100)
smpls <- unique(conditions)
cols <- c(1,rep(brewer.pal(8,"Dark2"),2))

### ==============================
## plot the first 3 dimensions
## The PCA looks good. first dimension
## is clearly developmental gradient
## that explains 37% of the data
### ==============================
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
              pch=c(rep(17,6),rep(19,12)))
legend("topleft",pch=c(17,17,19,19),col=pal,legend=smpls,bty="n")

scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
              pch=as.character(c(rep(1:3,2),rep(1:6,2))))
legend("topleft",fill=pal,legend=smpls,bty="n")

### ==============================
## Perform a Principal Component Analysis
## on the rld data
### ==============================
pc.rld <- prcomp(t(assay(rld)))
percent.rld <- round(summary(pc.rld)$importance[2,]*100)

### ==============================
## plot the first 3 dimensions
## The PCA looks good. first dimension
## is clearly developmental gradient
## that explains 37% of the data
### ==============================
scatterplot3d(pc.rld$x[,1],
              pc.rld$x[,2],
              pc.rld$x[,3],
              xlab=paste("Comp. 1 (",percent.rld[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent.rld[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent.rld[3],"%)",sep=""),
              color=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
              pch=c(rep(17,6),rep(19,12)))
legend("topleft",pch=c(17,17,19,19),col=pal,legend=smpls,bty="n")

### ==============================
## no big difference between RLT and VST
## going on with VST
### ==============================
par(mfrow=c(1,2))
scatterplot3d(pc.rld$x[,1],
              pc.rld$x[,2],
              pc.rld$x[,3],
              xlab=paste("Comp. 1 (",percent.rld[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent.rld[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent.rld[3],"%)",sep=""),
              color=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
              pch=as.character(c(rep(1:3,2),rep(1:6,2))))
legend("topleft",fill=c(0,pal),legend=c("RLT",smpls),bty="n")

scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
              pch=as.character(c(rep(1:3,2),rep(1:6,2))))
legend("topleft",fill=c(0,pal),legend=c("VST",smpls),bty="n")
par(mfrow=c(1,1))

### ==============================
## There seem to be more complex
## reason for the 2nd and 3rd
## dimensions. Clearly, the 
## difference between normal wood
## and tension wood, but as well
## the size of the developmental zones
## play a role there
## 1st - 2nd dimension
### ==============================
par(mar=c(5.1,4.1,4.1,2.1))
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
     pch=c(rep(17,6),rep(19,12)),
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft",pch=c(17,17,19,19),col=pal,legend=smpls,cex=0.8)

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
     pch=as.character(c(rep(1:3,2),rep(1:6,2))),
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft",fill=pal,legend=smpls,cex=0.8)

### ==============================
## 2nd -3rd dimension
### ==============================
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
     pch=c(rep(17,6),rep(19,12)),
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("topleft",pch=c(17,17,19,19),col=pal,legend=smpls,cex=0.8)

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=c(rep(pal[1:2],each=3),rep(pal[3:4],each=6)),
     pch=as.character(c(rep(1:3,2),rep(1:6,2))),
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("topleft",pch=c(17,17,19,19),col=pal,legend=smpls,cex=0.8)


### ==============================
## for K only
### ==============================

### ==============================
## For visualization, the data is
## submitted to a variance stabilization
## transformation
### ==============================
k.count.table <- count.table[,grep("k",colnames(count.table))]
conditions <- sub("-[1-6]","",colnames(k.count.table))
k.dds <- DESeqDataSetFromMatrix(countData = k.count.table,
                              colData = data.frame(condition=conditions),
                              design = ~ condition)
colData(k.dds)$condition <- factor(colData(k.dds)$condition,
                                 levels=unique(conditions))
k.vsd <- varianceStabilizingTransformation(k.dds, blind=TRUE)

### ==============================
## plot the mean agains the variance (sd)
## it should be flat if the VST worked properly
## it is not bu it is likely due to the 
## variability of the gene expression
## in normal wood compared to tension wood
## This is not a concern for visualization
## but might be for normalization
### ==============================
sum(rowSums(k.count.table)==0)/nrow(k.count.table)
meanSdPlot(assay(k.vsd)[rowSums(k.count.table)>0,], ylim = c(0,2.5))

value.list <- lapply(1:6,function(i,a){a[,i]},assay(k.vsd))
names(value.list) <- colnames(k.count.table)
plot.multidensity(value.list)

sel <- rowMeans(assay(k.vsd)) > 5
meanSdPlot(assay(k.vsd)[rowMeans(assay(vsd)) > 5,])

            

### ==============================
## Add some metrics in order to set
## cutoff for lowly-expressed genes
### ==============================
df <- as.data.frame(assay(vsd))
colnames(df) <- names(count.table)

## the vst adds a min value for non-expressed 
## gene. The following line of code is just used to 
## substract this value from all expression data.
## then all genes that are 0 can be considered 
## as non-expressed, which will be used for the
## analysis below 
df <- df - min(df)

ncols <- ncol(df)
df$expressed <- rowSums(df) > 0
df$mean.expression <- rowMeans(df[,1:ncols])
df$sample.expressed <- rowSums(df[,1:ncols]>0)
df$rmd <- apply(df[,1:ncols],1,rmd)

### ==============================
## Table of genes that are expressed 
## in all samples (TRUE) or not 
## expressed in all samples (FALSE)
## Here: 32685 genes are expressed
## in all samples
### ==============================

table(df$expressed)

### ==============================
## Table of genes that are expressed 
## in the given number of samples 
## (in none=0, in 1 sample=1 etc)
### ==============================
table(df$sample.expressed)

### ==============================
## Plot of mean expression accross 
## all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis)
### ==============================

#par(mar=c(3.1,2.1,2.1,0.1))
par(mar=c(4.1,4.1,1.1,0.1))
boxplot(split(df$mean.expression,df$sample.expressed))
title(ylab = "Mean expression", xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## Plot of mean variance accross 
## all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis)
### ==============================
boxplot(split(df$rmd,df$sample.expressed))
title(ylab = "Mean variance", xlab = "Number of samples with expressed genes", font.lab = 1)


### ==============================
## Plot the mean against the variance 
## after variance stabilization transformation
## and with mean expression values
## ordere from small to large
## X-axis indicates number of genes
## add two lines into the plot indicating
## 1. the first expressed gene 
## (in at least 1 sample) and
## 2. the first gene 
## expressed in all samples (this is the
## gene out of the 32685 all-sample expressed 
## ones with lowest mean expression
## accross all samples)
## Note that on the right side of this line
## there are still genes that are not 
## expressed in all samples (x-axis ordered)
## by mean expression, not by number of samples
## in which each gene is expressed
## We won't work with the ca. 18000 genes on
## on the left side of this line. But we will
## consider the 23000 genes on the right side
## of this line (! not necessarily expressed
## in all samples! )
### ==============================

meanSdPlot(assay(vsd))
nrow(assay(vsd))
abline(v=sum(!df$expressed))
abline(v=match("18",df[order(df$mean.expression),"sample.expressed"]))

### ==============================
## Get expression information for the
## gene that is indicated by the 
## right-most vertical line in the aforementioned
## plot. The mean expression value of
## this gene (1.424) will be our cut-off
## for considering genes that are (non-)expressed
## the rmd value here is Relative Mean Difference 
## (variance) and has a commons scale accross
## all genes
### ==============================

df[order(df$mean.expression),][match("18",df[order(df$mean.expression),"sample.expressed"]),]

### ==============================
## Now we set a threshold at that position
## and plot once again the mean expression 
## accross all samples (Y-axis) against
## number of samples in which
## these genes are expressed (x-axis) 
### ==============================

pos <- match("18",df[order(df$mean.expression),"sample.expressed"])
df <- df[order(df$mean.expression),]
dfs <- df[pos:nrow(df),]
boxplot(split(dfs$mean.expression,dfs$sample.expressed))
title(ylab = "Mean expression", xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## This plot below indicates
## relative mean difference (variance) 
## accross all samples against
## number of samples in which
## these genes are expressed (x-axis)
## Compared to what we had before 
## (when considering all genes), the
## variance has decreased and the 
## quality of the considered data
## therefore increased
### ==============================
boxplot(split(dfs$rmd,dfs$sample.expressed))
title(ylab = "Relative Mean Difference", xlab = "Number of samples with expressed genes", font.lab = 1)

### ==============================
## Looking more closely at the gene(s)
## that is/are expressed only in 5 or 6 
## of our samples we see that 
## the samples in which they are 
## expressed refer to similar and adjacent section
## positions in different biol. replicates
## This is very good! Shows that our
## strategy makes sense to filter out
## expression data that could be
## biologically relevant
### ==============================

dfs[grep("^5$",dfs$sample.expressed),]
dfs[grep("^6$",dfs$sample.expressed),]
dim(df)

## Save the file of normalized expression data for further analysis
## Don't forget to substract your threshold value (for us 1.424) before
## ongoing with data anlysis in excel. In our case we should end up with
## ca. 23000 expressed genes after that procedure. Minimum sample number
## with gene expression values should be 5 according to the last two plots.
write.csv(df,file="analysis/normalized-tensionwood-expression.csv")


#merge to annotation table there is a problem this doesn't work! GeneIDs not recognized in df2
annotationtable<-read.csv(file="/home/judith/ethylene-insensitive-HTSeq_DESeq/Ptrichocarpa_v3.0_210_gene-annotation_ERFs.csv")
annot_dataframe<-as.data.frame(annotationtable)
df <- read.csv(file="/home/judith/Manoj_TW_project_HTSeq/normalized-tensionwood-expression-RNASeqdata-dec2013.csv", stringsasfactors=FALSE)
rownames(df)<-as.character(df[,1])
cut.last<-function(text.in){
  temp<-NULL
  for(i in text.in){
    n.end<-nchar(i)
    temp<-c(temp,substr(i,1,(n.end-2)))
  }
  cut.last.out<-temp
  return(cut.last.out)
}
test<-cut.last(rownames(df))
rownames(df)<-test

df <- cbind(df,annot_dataframe[match(rownames(df),annot_dataframe$geneID),])
write.csv(df,"/home/judith/Manoj_TW_project_HTSeq/TWsectiondata.csv")
