#Trendy package
```{r}
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("Trendy")

library(Trendy)
```

#Data processing for TRENDY
```{r}
mRNAseqData=read.table("~/tmp_psychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
modsb1= mRNAseqData %>% 
  separate(Geneid,c("EnsemblID","Genename"),extra="merge")
modsb1$EnsemblID<-NULL
#Filtering for Akey alone
akeySestan1 <- modsb1 %>% filter(modsb1$Genename %in% results$hgnc_symbol) #From jm_Allen.R
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
```

#CEREBELLUM
```{r}
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

finalcbc <- t(cbcakey1)
colnames(finalcbc) <- time.vector 
set.seed(10)
res <- trendy(Data = finalcbc, tVectIn = time.vector, maxK = 3, minNumInSeg = 2, meanCut = 2)
res <- results(res)
res.top <- topTrendy(res)
#Good fit in Des & Pos sel.
rownames(res.top$Trends)[rownames(res.top$Trends) %in% both$hgnc_symbol]

res.trend <- trendHeatmap(res.top)
##Heatmap homemad
class(res.top$Trends)
pheatmap(res.top$Trends, cluster_rows=FALSE, cluster_cols = FALSE, legend_breaks = c(-1,0,1), filename = "~/raul_tesina/2.plots/trendy_segmentedRegression/CBC_UpDown_Deserts.pdf")

sb_inDesPey <- res.top$Trends[rownames(res.top$Trends)[rownames(res.top$Trends) %in% both$hgnc_symbol],]

pheatmap(sb_inDesPey, cluster_rows=FALSE, cluster_cols = FALSE, legend_breaks = c(-1,0,1),filename = "~/raul_tesina/2.plots/trendy_segmentedRegression/CBC_UpDown_DesertsPosSel.pdf")


#Selected genes
par(mfrow=c(3,2))
plot2 <- plotFeature(finalcbc,tVectIn = time.vector,featureNames = c("SYT6", "ROBO2", "CADPS2", "GPR22", "BCAP29"),trendyOutData = res)

pdf("~/raul_tesina/2.plots/trendy_segmentedRegression/CBC_GenesDesertPosSel.pdf")
par(mfrow=c(3,2))
plot2 <- plotFeature(finalcbc,tVectIn = time.vector,featureNames = c("SYT6", "ROBO2", "CADPS2", "GPR22", "BCAP29"),trendyOutData = res)
dev.off()
```

#STRIATUM
```{r}
strakey <- filter(finalakeySestan1, Regioncode == "STR")
strakey <- arrange(strakey, Window)
time.vector <- strakey$Window
strakey <- strakey %>% unite("code", c("Braincode","Regioncode","Window"))
rownames(strakey) <- strakey$code
strakey$code <- NULL

strakey1 <- as_tibble(strakey)

cols.num <- colnames(strakey1)
strakey1[,cols.num] <- lapply(strakey1[cols.num],as.character)
strakey1[,cols.num] <- lapply(strakey1[cols.num],as.numeric)

finalstr <- t(strakey1)
colnames(finalstr) <- time.vector 
set.seed(10)
res <- trendy(Data = finalstr, tVectIn = time.vector, maxK = 3, minNumInSeg = 2, meanCut = 2)
res <- results(res)
res.top <- topTrendy(res)
#Good fit in Des & Pos sel.
rownames(res.top$Trends)[rownames(res.top$Trends) %in% both$hgnc_symbol]

res.trend <- trendHeatmap(res.top)
##Heatmap homemad
class(res.top$Trends)
pheatmap(res.top$Trends, cluster_rows=FALSE, cluster_cols = FALSE, legend_breaks = c(-1,0,1), filename = "~/raul_tesina/2.plots/trendy_segmentedRegression/STR_UpDown_Deserts.pdf")

sb_inDesPey <- res.top$Trends[rownames(res.top$Trends)[rownames(res.top$Trends) %in% both$hgnc_symbol],]

pheatmap(sb_inDesPey, cluster_rows=FALSE, cluster_cols = FALSE, legend_breaks = c(-1,0,1),filename = "~/raul_tesina/2.plots/trendy_segmentedRegression/STR_UpDown_DesertsPosSel.pdf")


#Selected genes
par(mfrow=c(3,2))
plot2 <- plotFeature(finalstr,tVectIn = time.vector,featureNames = c("ROBO2", "BCAP29", "ST7"),trendyOutData = res)

pdf("~/raul_tesina/2.plots/trendy_segmentedRegression/STR_GenesDesertPosSel.pdf")
par(mfrow=c(3,2))
plot2 <- plotFeature(finalstr,tVectIn = time.vector,featureNames = c("ROBO2", "BCAP29", "ST7"),trendyOutData = res)
dev.off()
```

#MD
```{r}
mdakey <- filter(finalakeySestan1, Regioncode == "MD")
mdakey <- arrange(mdakey, Window)
time.vector <- mdakey$Window
mdakey <- mdakey %>% unite("code", c("Braincode","Regioncode","Window"))
rownames(mdakey) <- mdakey$code
mdakey$code <- NULL

mdakey1 <- as_tibble(mdakey)

cols.num <- colnames(mdakey1)
mdakey1[,cols.num] <- lapply(mdakey1[cols.num],as.character)
mdakey1[,cols.num] <- lapply(mdakey1[cols.num],as.numeric)

finalmd <- t(mdakey1)
colnames(finalmd) <- time.vector 
set.seed(10)
res <- trendy(Data = finalmd, tVectIn = time.vector, maxK = 3, minNumInSeg = 2, meanCut = 2)
res <- results(res)
res.top <- topTrendy(res)
#Good fit in Des & Pos sel.
rownames(res.top$Trends)[rownames(res.top$Trends) %in% both$hgnc_symbol]

res.trend <- trendHeatmap(res.top)
##Heatmap homemad
class(res.top$Trends)
pheatmap(res.top$Trends, cluster_rows=FALSE, cluster_cols = FALSE, legend_breaks = c(-1,0,1), filename = "~/raul_tesina/2.plots/trendy_segmentedRegression/MD_UpDown_Deserts.pdf")

sb_inDesPey <- res.top$Trends[rownames(res.top$Trends)[rownames(res.top$Trends) %in% both$hgnc_symbol],]

pheatmap(sb_inDesPey, cluster_rows=FALSE, cluster_cols = FALSE, legend_breaks = c(-1,0,1),filename = "~/raul_tesina/2.plots/trendy_segmentedRegression/MD_UpDown_DesertsPosSel.pdf")


#Selected genes
par(mfrow=c(3,2))
plot2 <- plotFeature(finalmd,tVectIn = time.vector,featureNames = c("ROBO2", "BCAP29", "ST7"),trendyOutData = res)

pdf("~/raul_tesina/2.plots/trendy_segmentedRegression/MD_GenesDesertPosSel.pdf")
par(mfrow=c(3,2))
plot2 <- plotFeature(finalmd,tVectIn = time.vector,featureNames = c("ROBO2", "BCAP29", "ST7"),trendyOutData = res)
dev.off()
```
