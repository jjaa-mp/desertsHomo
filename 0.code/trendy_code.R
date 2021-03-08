#Trendy package
```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("Trendy")

library(Trendy)
```

#Example with deserts of introgression and cerebellum
```{r}
mRNAseqData=read.table("~/tmp_psychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
modsb1= mRNAseqData %>% 
  separate(Geneid,c("EnsemblID","Genename"),extra="merge")
modsb1$EnsemblID<-NULL
#Filtering for Akey alone
akeySestan1 <- modsb1 %>% filter(modsb1$Genename %in% results$hgnc_symbol)
akeySestan1=t(akeySestan1)

#As dataframe
akeySestan1=as.data.frame(akeySestan1)

colnames(akeySestan1) <- as.matrix(unlist(akeySestan1[1,]))
akeySestan1 <- akeySestan1[-1, ]

akeySestan1 <- cbind(info = rownames(akeySestan1), akeySestan1)
rownames(akeySestan1) <- 1:nrow(akeySestan1)

#duplicated columns - issue in raw data. Here: 
colnames(akeySestan1)[duplicated(colnames(akeySestan1))] #1
akeySestan1 <- akeySestan1[, !duplicated(colnames(akeySestan1))]

akeySestan1=akeySestan1 %>% separate(info, c("Braincode","Regioncode"))


#Brining the metadata of the database
library(readxl)
metadatamRNAseq=read_xlsx("~/tmp_psychENCODE/mRNA-seq_QC.xlsx",skip = 3)
modMetadatamRNAseq=na.omit(metadatamRNAseq)
modMetadatamRNAseq=modMetadatamRNAseq %>% select(1:3)
modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)

finalakeySestan1=merge(modMetadatamRNAseq,akeySestan1,by=c("Braincode", "Regioncode"))

cbcakey <- filter(finalakeySestan1, Regioncode == "CBC")

cbcakey$Braincode <- NULL
cbcakey$Regioncode <- NULL

cbcakey1 <- as_tibble(cbcakey)

cols.num <- colnames(cbcakey1)
cbcakey1[,cols.num] <- lapply(cbcakey1[cols.num],as.character)
cbcakey1[,cols.num] <- lapply(cbcakey1[cols.num],as.numeric)

cbcakey1_summ <- cbcakey1 %>%
  group_by(Window) %>%
  dplyr::summarise_all(median, na.rm=TRUE)

finalcbc <- t(cbcakey1_summ)
nms <- finalcbc[1,]
colnames(finalcbc) <- nms

finalcbc<-finalcbc[-1,] #READY FOR TRENDY


time.vector <- 2:9
res <- trendy(Data = finalcbc, tVectIn = time.vector, maxK = 3, minNumInSeg = 2)
res <- results(res)
res.top <- topTrendy(res)
res.trend <- trendHeatmap(res.top)
library(gplots)
heatmap.2(finalcbc[names(res.trend$firstup),],trace="none", Rowv=FALSE,Colv=FALSE,dendrogram='none',scale="row", main="top genes (first go up)")

par(mfrow=c(3,2))
plotFeature(Data = finalcbc, tVectIn = time.vector, simple = TRUE,featureNames = names(res.trend$firstnochange)[1],trendyOutData = res)

par(mfrow=c(3,2))#specify the layout of multiple plots in a single panel
plotFeature(Data = finalcbc, tVectIn = time.vector, simple = FALSE,showLegend = TRUE, legendLocation='side',cexLegend=1,featureNames = names(res.trend$firstnochange)[13],trendyOutData = res)
```

#Alternative keeping replicates
```{r}
mRNAseqData=read.table("~/tmp_psychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
modsb1= mRNAseqData %>% 
  separate(Geneid,c("EnsemblID","Genename"),extra="merge")
modsb1$EnsemblID<-NULL
#Filtering for Akey alone
akeySestan1 <- modsb1 %>% filter(modsb1$Genename %in% results$hgnc_symbol)
akeySestan1=t(akeySestan1)

#As dataframe
akeySestan1=as.data.frame(akeySestan1)

colnames(akeySestan1) <- as.matrix(unlist(akeySestan1[1,]))
akeySestan1 <- akeySestan1[-1, ]

akeySestan1 <- cbind(info = rownames(akeySestan1), akeySestan1)
rownames(akeySestan1) <- 1:nrow(akeySestan1)

#duplicated columns - issue in raw data. Here: 
colnames(akeySestan1)[duplicated(colnames(akeySestan1))] #1
akeySestan1 <- akeySestan1[, !duplicated(colnames(akeySestan1))]

akeySestan1=akeySestan1 %>% separate(info, c("Braincode","Regioncode"))


#Brining the metadata of the database
library(readxl)
metadatamRNAseq=read_xlsx("~/tmp_psychENCODE/mRNA-seq_QC.xlsx",skip = 3)
modMetadatamRNAseq=na.omit(metadatamRNAseq)
modMetadatamRNAseq=modMetadatamRNAseq %>% select(1:3)
modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)

finalakeySestan1=merge(modMetadatamRNAseq,akeySestan1,by=c("Braincode", "Regioncode"))

cbcakey <- filter(finalakeySestan1, Regioncode == "CBC")
cbcakey <- arrange(cbcakey, Window)
time.vector <- cbcakey$Window
cbcakey <- cbcakey %>% unite("code", c("Braincode","Regioncode","Window"))
rownames(cbcakey) <- cbcakey$code
cbcakey$code <- NULL

cbcakey1 <- as_tibble(cbcakey)

cols.num <- colnames(cbcakey1)
cbcakey1[,cols.num] <- lapply(cbcakey1[cols.num],as.character)
cbcakey1[,cols.num] <- lapply(cbcakey1[cols.num],as.numeric)

#cbcakey1_summ <- cbcakey1 %>%
  #group_by(Window) %>%
  #dplyr::summarise_all(median, na.rm=TRUE)

finalcbc <- t(cbcakey1)
colnames(finalcbc) <- time.vector #READY FOR TRENDY
saveRDS(finalcbc, file = "~/CBC_trendy.rds")

set.seed(10)
res <- trendy(Data = finalcbc, tVectIn = time.vector, maxK = 3, minNumInSeg = 2, meanCut = 2)
res <- results(res)
res.top <- topTrendy(res)
#Good fit in Des & Pos sel.
rownames(res.top$Trends)[rownames(res.top$Trends) %in% both$hgnc_symbol]

res.trend <- trendHeatmap(res.top)
##Heatmap homemad
class(res.top$Trends)
pheatmap(res.top$Trends, cluster_rows=FALSE, cluster_cols = FALSE, legend_breaks = c(-1,0,1), filename = "~/raul_tesina/2.plots/trendy_segmentedRegression/pheatmap_UpDown_Deserts.pdf")

sb_inDesPey <- res.top$Trends[rownames(res.top$Trends)[rownames(res.top$Trends) %in% both$hgnc_symbol],]
pheatmap(sb_inDesPey, cluster_rows=FALSE, cluster_cols = FALSE, legend_breaks = c(-1,0,1),filename = "~/raul_tesina/2.plots/trendy_segmentedRegression/pheatmap_UpDown_DesertsPosSel.pdf")

#library(gplots)
#heatmap.2(finalcbc[names(res.trend$firstup),],trace="none", Rowv=FALSE,Colv=FALSE,dendrogram='none',scale="row", main="top genes (first go up)")


#Individual genes:
names(res.trend$firstnochange)
names(res.trend$firstnochange)[names(res.trend$firstnochange) %in% both$hgnc_symbol]

par(mfrow=c(3,2))
plotFeature(Data = finalcbc, tVectIn = time.vector, simple = TRUE,featureNames =names(res.trend$firstnochange)[c(1,6,8)],trendyOutData = res)

par(mfrow=c(3,2))
plotFeature(Data = finalcbc, tVectIn = time.vector, simple = FALSE,showLegend = FALSE, featureNames = names(res.trend$firstnochange)[c(1,6,8)],trendyOutData = res)


#Genes with specific patterns
##Genes that increase after first stages (e.g CADPS2)
extractPattern(res, Pattern = c("same", "up"))
par(mfrow=c(3,2))
plotFeature(Data = finalcbc, tVectIn = time.vector, simple = FALSE,showLegend = FALSE, featureNames = names(res.trend$firstnochange)[c(6,8,9,13)],trendyOutData = res)

#Selected genes
par(mfrow=c(3,2))
plot2 <- plotFeature(finalcbc,tVectIn = time.vector,featureNames = "CELSR2",trendyOutData = res)

pdf("~/raul_tesina/2.plots/trendy_segmentedRegression/trendy_GenesDesertPosSel.pdf")
par(mfrow=c(3,2))
plot2 <- plotFeature(finalcbc,tVectIn = time.vector,featureNames = c("SYT6", "ROBO2", "CADPS2", "GPR22", "BCAP29"),trendyOutData = res)
dev.off()

#Timepoints breaks - all genes
res.bp <- breakpointDist(res.top)
barplot(res.bp, ylab="Number of breakpoints", col="lightblue")
```
