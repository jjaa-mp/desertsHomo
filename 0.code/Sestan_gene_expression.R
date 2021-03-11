#Packages
library(biomaRt)
library(data.table)
library(ggplot2)
library(dplyr)
library(GGally)
library(viridis)
library(ggpubr)
library(grid)
library(gridExtra)
library(lattice)
library(xlsx)
library(tidyr)
library(tidyverse)
library(DescTools)
library(reshape2)
library(dplyr)
library(Matrix)
library(pheatmap)
library(plyr)
library(rstatix)


# mRNAseqData=read.table("mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
# modsb1= mRNAseqData %>% 
#   separate(Geneid,c("EnsemblID","Genename"),extra="merge")
# modsb1$EnsemblID<-NULL


crossedDataRNAseq<-read.csv("crossedDataBrainPsychmRNAseq.csv")
crossedDataRNAseq$X<-NULL

#Transform to a usefull matrix

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- row.names(DataRNAseq[[2]])
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                                                          "hgnc_symbol", "description"),values=genes,mart= mart)


row.names(DataRNAseqw2Fin)
DataRNAseq<-list()
#windows
for(i in 2:9){
  DataRNAseqw<-crossedDataRNAseq %>% filter(crossedDataRNAseq$Window==i)
  DataRNAseqw$FullReg<-NULL
  DataRNAseqw<-DataRNAseqw[,8:length(DataRNAseqw)]
  DataRNAseqw<-t(DataRNAseqw)
  colnames(DataRNAseqw) <- DataRNAseqw[1,]
  DataRNAseqw <- DataRNAseqw[-1, ]
  DataRNAseqw<-as.data.frame(DataRNAseqw)
  nms<-unique(names(DataRNAseqw))
  DataRNAseqw[] <- as.data.frame(lapply(DataRNAseqw, function(x) as.numeric(as.character(x))))
  
  DataRNAseq[[i]]<-as.data.frame(sapply(nms, function(x)  rowMeans(as.data.frame(DataRNAseqw)[names(as.data.frame(DataRNAseqw)) %in% x])))
  # DataRNAseq[["2"]]<-t(rowsum(t(DataRNAseqw), names(DataRNAseqw))/c(table(names(DataRNAseqw))))
}

#hgnc_symbol
for(i in 2:9){
  DataRNAseq[[i]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(DataRNAseq[[i]]), G_list$ensembl_gene_id)]
}


#Filtering for Akey alone
akeySestan1 <- modsb1 %>% filter(modsb1$Genename %in% results$hgnc_symbol)
# akeySestan1=t(akeySestan1)
cpaleySestan<-data.frame(akeySestan1)#copy
cpaleySestan$Genename<-NULL
meanAkeySestan<-data.frame(genename=akeySestan1$Genename,meanexp=rowMeans(cpaleySestan))

q75Gen<-c()
q75genes<-meanAkeySestan %>% filter(meanAkeySestan$meanexp > quantile(meanAkeySestan$meanexp, 0.75))




for (h in 2:length(names(DataRNAseq[[i]]))){
  q75 <- ab[[i]][h] %>% filter(ab[[i]][h] > quantile(ab[[i]][[h]], 0.75))
}


DataRNAseqTemp<-list()
q75 = vector(mode="list", length = length(DataRNAseq))
for (i in 2:length(DataRNAseq)){
  DataRNAseqTemp[[i]] <- DataRNAseq[[i]] %>% filter(rownames(DataRNAseq[[i]]) %in% results$ensembl_gene_id)
  for (h in 1:(length(names(DataRNAseqTemp[[i]]))-1)){
    #Akey
    q75Temp<-as.data.frame(DataRNAseqTemp[[i]] %>% arrange("gene_name",names(DataRNAseqTemp[[i]])[h]) %>% filter(DataRNAseqTemp[[i]][h] > quantile(DataRNAseqTemp[[i]][[h]], 0.75)))[,c("gene_name",names(DataRNAseqTemp[[i]])[h])]
    q75[[i]][[names(DataRNAseqTemp[[i]])[h]]] <- q75Temp
  }
}
install.packages("xlsx")
require(openxlsx)
for (i in 2:length(DataRNAseq)){
  for (h in 1:length(q75[[i]])){
    if (dim(q75[[i]][[names(q75[[i]])[h]]])[1] == 0) next #To skip empty dataframes
    write.xlsx(q75[[i]], file=paste("Sestan_Akey_highExprq75_w",toString(i),".xlsx"),append=TRUE) #Add age number at the beginning of the sheet name
  }
}
