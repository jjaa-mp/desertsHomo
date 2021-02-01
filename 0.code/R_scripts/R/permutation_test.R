permutation_test <- function(npermutations, abdesert) {
  #This function takes abadult & the deserts coordinates from biominput() as input
  
  mask <- data.frame(c("chr1","chr3","chr7","chr8"), 
                     c(105400000, 74100000,106200000,49400000), 
                     c(120600000, 89300000,123200000,66500000))
  
  #More efficient this way
  random <- replicate(n = npermutations, get_data_perms(mask), simplify = FALSE)
  test <- lapply(random, randomregion_biomart_query)
  
  #get the expression data of random genes
  loadd(abadult)
  #Internal function
  selected <- lapply(test, filter_genexpr_random)
  # quantile 10 cutoff to avoid 0's happens here

  # Initializing calculation of means (of log) per substructures 
  randomregions_struc <- lapply(selected, means_log)
  # log transformation happens hear
  
  #randomregions_struc now contains:
  #the means of gene expresion (normalized with log, cutoff included)
  # per structure, per random region generated
  # so 1000 permutations = 1000 nested tables
  
  #Now,
  loadd(abakey_data)
} 
