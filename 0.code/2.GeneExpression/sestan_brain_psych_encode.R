library("readxl")
library("dplyr")
library("tidyr")
mRNAseqData=read.table("mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)

#Keeping just one gene identifier
modmRNAseqData=mRNAseqData %>% 
  separate(Geneid,c("Geneid","NormalName"),extra="merge")
modmRNAseqData$NormalName<-NULL
drop(modmRNAseqData$NormalName)
#Transpose the matrix to cross it with other data
modmRNAseqData=t(modmRNAseqData)

#Transformation of the original data

modmRNAseqData=as.data.frame(modmRNAseqData)

colnames(modmRNAseqData) <- as.matrix(unlist(modmRNAseqData[1,]))
modmRNAseqData <- modmRNAseqData[-1, ]

modmRNAseqData <- cbind(info = rownames(modmRNAseqData), modmRNAseqData)
rownames(modmRNAseqData) <- 1:nrow(modmRNAseqData)

modmRNAseqData=modmRNAseqData %>% 
  separate(info,c("Braincode","regioncode"))
staticmodmRNAseqData=modmRNAseqData

# Reducing the data to work with it for the crossdata

redmodmRNAseqData=head(modmRNAseqData)
redmodmRNAseqData=staticredmodmRNAseqData
colnames(redmodmRNAseqData)[1]<-"Braincode"

#Brining the metadata of the database
metadatamRNAseq=read_xlsx("mRNA-seq_Sample metadata.xlsx",skip = 3)


modMetadatamRNAseq=na.omit(metadatamRNAseq)
modMetadatamRNAseq=modMetadatamRNAseq %>% select(1:7)
modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)

#Importing the brain region acronims correspondance
brainRegionCorrespond=read.table("brainRegionCorresp.csv",sep=":",header = TRUE)
prmodmRNAseqData=merge(brainRegionCorrespond,modmRNAseqData,by="regioncode")
# prredmodmRNAseqData=merge(brainRegionCorrespond,redmodMetadatamRNAseq,by="regioncode")

finalredmRNAseqData=merge(modMetadatamRNAseq,prmodmRNAseqData,by="Braincode")
# testredmRNAseqData=merge(modMetadatamRNAseq,prredmodmRNAseqData,by="Braincode")

write.csv(finalredmRNAseqData,"crossedDatamRNAseq.csv")

