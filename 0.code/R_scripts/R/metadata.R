
metadata <- function(logsestan){
  #Brining the metadata of the database
  # This path is specific of my computer: should have download of dataframe integrated
  # (plus it breaks github - space limit)
  metadatamRNAseq=read_xlsx("../../../rawsestan/mRNA-seq_QC.xlsx",skip = 3)
  
  
  modMetadatamRNAseq=na.omit(metadatamRNAseq)
  modMetadatamRNAseq=modMetadatamRNAseq %>% dplyr::select(1:3)
  modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)
  
  #With log
  finalsestan=merge(modMetadatamRNAseq,logsestan,by=c("Braincode", "Regioncode"))
  #Without log:
  #finalsestan=merge(modMetadatamRNAseq,akeypeySestan,by=c("Braincode", "Regioncode"))
  
  finalsestan$Braincode <- NULL
  finalsestan$Label<- paste(finalsestan$Regioncode, finalsestan$Window, sep="_")
  finalsestan$Regioncode <- NULL
  finalsestan$Window <- NULL
  
  finalsestan <- as_tibble(finalsestan)
  
  finalsestan_w_mean <- finalsestan %>%
    group_by(Label) %>%
    dplyr::summarise_all(median, na.rm=TRUE)
  
  #For plotting
  dfsestan <- data.frame(Structure=finalsestan_w_mean$Label, Means=rowMeans(finalsestan_w_mean[,-1]))
  
  dfsestan <- dfsestan %>% 
    separate(Structure, c("Structure","Window"))
  preli <- pivot_wider(dfsestan, names_from = Window, values_from = Means)
  preli[,10] <- NULL
  preli <- preli[complete.cases(preli), ]
  colnames(preli) <- c("Structure", "Fetal1", "Fetal2", "Fetal3", "Birth/Infant", "Infan/Childh", "Childh", "Adolescence", "Adulth")
  #PLOT
  levels(colnames(preli)) <- c("Structure", "Fetal1", "Fetal2", "Fetal3", "Birth/Infant", "Infan/Childh", "Childh", "Adolescence", "Adulth")
  return(preli)
}