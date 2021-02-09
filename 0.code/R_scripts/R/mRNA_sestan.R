 mRNA_sestan <- function(region) {
  # This path is specific of my computer: should have download of data integrated
  # (plus it breaks github - space limit)
  mRNAseqData=read.table("../../../rawsestan/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
  modsb1= mRNAseqData %>% 
    separate(Geneid,c("EnsemblID","Genename"),extra="merge")
  modsb1$EnsemblID<-NULL
  #Filtering for Akey alone

    sestan <- modsb1 %>% filter(modsb1$Genename %in% region$hgnc_symbol)
    sestan <- t(sestan)

  
  #As dataframe
  sestan <- as.data.frame(sestan)
  
  colnames(sestan) <- as.matrix(unlist(sestan[1,]))
  sestan <- sestan[-1, ]
  
  sestan <- cbind(info = rownames(sestan), sestan)
  rownames(sestan) <- 1:nrow(sestan)
  
  #duplicated columns - issue in raw data. Here: 
  colnames(sestan)[duplicated(colnames(sestan))] #1
  sestan <- sestan[, !duplicated(colnames(sestan))]
  
  sestan=sestan %>% 
    separate(info, c("Braincode","Regioncode"))
  
  #Normality
  #normcheck <- akeypeySestan %>% pivot_longer(!("Braincode"|"Regioncode"), names_to = "Genes", values_to = "RPKM")
  #normcheck$RPKM <- as.numeric(normcheck$RPKM)
  #ggqqplot(normcheck, "RPKM", facet.by = "Regioncode")
  
  #Transformation
  logsestan <- sestan
  
  cols.num <- colnames(logsestan[3:ncol(logsestan)])
  
  logsestan[,cols.num] <- lapply(logsestan[cols.num],as.character)
  logsestan[,cols.num] <- lapply(logsestan[cols.num],as.numeric)
  
  logsestan[3:ncol(logsestan)] <- log2(logsestan[3:ncol(logsestan)]+1)
  
  logsestan <- metadata(logsestan)
  return(logsestan)
 }
 
  #Normality after transformation
  #lognormcheck <- logakeypeySestan %>% pivot_longer(!("Braincode"|"Regioncode"), names_to = "Genes", values_to = "RPKM")
  #lognormcheck$RPKM <- as.numeric(lognormcheck$RPKM)
  #ggqqplot(lognormcheck, "RPKM", facet.by = "Regioncode")

