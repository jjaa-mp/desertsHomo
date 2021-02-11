clean_sestan <- function(){
  cleaningdf <- read.csv("../../../rawsestan/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",
                         sep = "\t")
  #Needs an special function because of data structure
  cleaningdf <- cleaningdf %>% 
    separate(Geneid,c("EnsemblID","gene_name"),extra="merge")
  cleaningdf$EnsemblID<-NULL
  return(cleaningdf)
} 
