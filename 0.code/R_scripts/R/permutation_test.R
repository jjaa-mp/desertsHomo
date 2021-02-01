permutation_test <- function(abtemp, akey) {
  #This function takes abadult & the deserts coordinates from biominput() as input
  
  mask <- data.frame(c("chr1","chr3","chr7","chr8"), 
                     c(105400000, 74100000,106200000,49400000), 
                     c(120600000, 89300000,123200000,66500000))
  
  #More efficient this way
  random <- replicate(n = 10, get_data_perms(mask), simplify = FALSE)
  test <- lapply(random, randomregion_biomart_query)
  
  # Check if the genes AB data are in deserts of introgression
  abDeserts <- abtemp[abtemp$gene_name %in% akey$hgnc_symbol,]
  
  # drop gene IDs
  abDeserts <- NULL
  drop(abDeserts$gene_name)
  
  # Initializing calculation of means per substructures with the proper dataframe
  meansSubstruct <- colMeans(abDeserts)
  rownames(meansSubstruct)
  meansSubstructTemp <- cbind(idBrainStruct = rownames(meansSubstruct), meansSubstruct)
  meansSubstructTemp <- cbind(idBrainStruct = rownames(meansSubstructTemp), meansSubstructTemp)
  rownames(meansSubstructTemp) <- 1:nrow(meansSubstructTemp)
  
  
  #abtempFin<-abtemp2[abtemp2$idGene %in% resultsTemp$ensembl_gene_id,]
  
  
} 
