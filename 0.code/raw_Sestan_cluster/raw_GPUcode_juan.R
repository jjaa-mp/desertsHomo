#!usr/bin/env Rscript

library(magrittr)
library(dplyr)
library(tidyr)
library(readxl)

mRNAseqData <- read.table("mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)

modsb1<- mRNAseqData %>% 
  separate(Geneid,c("EnsemblID","Genename"),extra="merge")
modsb1$EnsemblID<-NULL
#Filtering for Akey alone
rawSestan <- modsb1
rm(modsb1)

rawSestan <- t(rawSestan)

#As dataframe
rawSestan <- as.data.frame(rawSestan)

colnames(rawSestan) <- as.matrix(unlist(rawSestan[1,]))
rawSestan <- rawSestan[-1, ]

rawSestan <- cbind(info = rownames(rawSestan), rawSestan)
rownames(rawSestan) <- 1:nrow(rawSestan)

#duplicated columns - Remove duplicates if needed
colnames(rawSestan)[duplicated(colnames(rawSestan))] #0
rawSestan <- rawSestan[, !duplicated(colnames(rawSestan))]

rawSestan <- rawSestan %>% 
  separate(info, c("Braincode","Regioncode"))

# #Normality
# normcheck <- rawSestan %>% pivot_longer(!("Braincode"|"Regioncode"), names_to = "Genes", values_to = "RPKM")
# normcheck$RPKM <- as.character(normcheck$RPKM)
# normcheck$RPKM <- as.numeric(normcheck$RPKM)
# 
# norm_1=merge(modMetadatamRNAseq,normcheck,by=c("Braincode", "Regioncode"))
# 
# norm_1 <- as_tibble(norm_1)
# 
# ggqqplot(norm_1, "RPKM", facet.by = "Window")

 
# #Transformation
lograwsestan1 <- rawSestan
cols.num <- colnames(lograwsestan1[3:ncol(lograwsestan1)])
lograwsestan1[,cols.num] <- lapply(lograwsestan1[cols.num],as.character)
lograwsestan1[,cols.num] <- lapply(lograwsestan1[cols.num],as.numeric)
lograwsestan1[3:ncol(lograwsestan1)] <- log2(lograwsestan1[3:ncol(lograwsestan1)]+1)

#Brining the metadata of the database
metadatamRNAseq <- read_xlsx("mRNA-seq_QC.xlsx",skip = 3)

modMetadatamRNAseq <- na.omit(metadatamRNAseq)
modMetadatamRNAseq <- modMetadatamRNAseq %>% select(1:3)
modMetadatamRNAseq <- as.data.frame(modMetadatamRNAseq)

#With log
finalrawSestan <- merge(modMetadatamRNAseq,lograwsestan1,by=c("Braincode", "Regioncode"))
#Without log:
#finalrawSestan=merge(modMetadatamRNAseq,rawSestan,by=c("Braincode", "Regioncode"))
rm(modMetadatamRNAseq, lograwsestan1)

finalrawSestan$Braincode <- NULL
finalrawSestan$Label<- paste(finalrawSestan$Regioncode, finalrawSestan$Window, sep="_")
finalrawSestan$Regioncode <- NULL
finalrawSestan$Window <- NULL

finalrawSestan <- as_tibble(finalrawSestan)

finalrawSestan_w_mean <- finalrawSestan %>%
    group_by(Label) %>%
    dplyr::summarise_all(median, na.rm=TRUE)


#For plotting
dfrawSestan <- data.frame(Structure=finalrawSestan_w_mean$Label, Means=rowMeans(finalrawSestan_w_mean[,-1]))
rm(finalrawSestan_w_mean, finalrawSestan)

dfrawSestan <- dfrawSestan %>% 
  separate(Structure, c("Structure","Window"))


df_raw <- pivot_wider(dfrawSestan, names_from = Window, values_from = Means)
df_raw[,10] <- NULL
df_raw <- df_raw[complete.cases(df_raw), ]
colnames(df_raw) <- c("Structure", "Fetal1", "Fetal2", "Fetal3", "Birth/Ifan", "Infan/Childh", "Childh", "Adolescence", "Adulth")
#PLOT
levels(colnames(df_raw)) <- c("Structure", "Fetal1", "Fetal2", "Fetal3", "Birth/Ifancy", "Infancy/Childh", "Childh", "Adolescence", "Adulth")


write.csv(df_raw, "df_raw.csv")
