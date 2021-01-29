get_data_perms <- function(abtemp, akey) { 
  #takes abadult, deserts corrdinates from biominput()
  
  #Extracted gene names from coordinate regions via bioMart in previous step

  abtemp 
  #same with adult data BA
  
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
  
  
  # Permutations
  # create mask:
  mask <- data.frame(c("chr1","chr3","chr7","chr8"), 
                     c(105400000, 74100000,106200000,49400000), 
                     c(120600000, 89300000,123200000,66500000))
  randReg50set<-createRandomRegions(nregions = 50,length.mean = 15000000,
                                     length.sd = 1000000,genome = "hg19", mask = mask)

  randomdesert<-c(runValue(seqnames(randReg50set)),start(randReg50set),end(randReg50set))
  # See comments in permutations code as to why this is implemented this way 
  #work in progress
  
  abtempFin<-abtemp2[abtemp2$idGene %in% resultsTemp$ensembl_gene_id,]

}