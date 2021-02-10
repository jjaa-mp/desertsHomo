permutation_test_sestan <- function(input) {
  
  cleaningdf <- read.csv("../../../rawsestan/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",
                    sep = "\t")
  #Needs an special function because of data structure
  cleaningdf <- cleaningdf %>% 
    separate(Geneid,c("EnsemblID","Genename"),extra="merge")
  cleaningdf$EnsemblID<-NULL
  #Filtering for Akey alone

  cleaningdf <- t(cleaningdf)
  
  cleaningdf <- as.data.frame(cleaningdf)
  cleaningdf <- tibble::rownames_to_column(cleaningdf, "structures")
  
  colnames(cleaningdf) <- cleaningdf[1,]
  
  cleaningdf <- cleaningdf[-1,]
  
  #duplicated columns - issue in raw data. Here: 
  colnames(cleaningdf)[duplicated(colnames(cleaningdf))] #1
  cleaningdf <- cleaningdf[, !duplicated(colnames(cleaningdf))]
  
  cleaningdf <- cleaningdf %>% 
    separate(Genename, c("Braincode","Regioncode"))
  
  metadatamRNAseq <- read_xlsx("../../../rawsestan/mRNA-seq_QC.xlsx",skip = 3)
  
  
  modMetadatamRNAseq <- na.omit(metadatamRNAseq)
  modMetadatamRNAseq <- modMetadatamRNAseq %>% dplyr::select(1:3)
  modMetadatamRNAseq <- as.data.frame(modMetadatamRNAseq)
  
  #With log
  cleaningdf <- merge(modMetadatamRNAseq,cleaningdf,by=c("Braincode", "Regioncode"))
  
  cleaningdf$Braincode <- NULL
  cleaningdf$Label<- paste(as.character(cleaningdf$Regioncode), as.character(cleaningdf$Window), sep="_")
  cleaningdf$Regioncode <- NULL
  cleaningdf$Window <- NULL
  
  cleaningdf <- as_tibble(cleaningdf)
  cleaningdf$Label <- as.factor(cleaningdf$Label)

  cleaningdf <- cleaningdf %>%
    group_by(as.character(cleaningdf$Label)) %>%
    dplyr::summarise_all(median, na.rm=TRUE)
  
}