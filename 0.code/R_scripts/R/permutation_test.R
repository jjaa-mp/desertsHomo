permutation_test <- function(npermutations, dataset) {
  #This function takes abadult & the deserts coordinates from biominput() as input
  # Outputs random regions means by structure to compare against expression in special subsets
  # won't work if lock_envir = TRUE in _drake.R
  
  mask <- data.frame(c("chr1","chr3","chr7","chr8"), 
                     c(105400000, 74100000,106200000,49400000), 
                     c(120600000, 89300000,123200000,66500000))
  
  #More efficient this way
  random <- replicate(n = npermutations, get_data_perms(mask), simplify = FALSE)
  run <- lapply(random, randomregion_biomart_query)
  #remove empty dataframes
  cleanrun <- run %>% 
    purrr::discard(empty)
  
  #Internal functions
  selected <- lapply(cleanrun, filter_genexpr_random, abadult = dataset)
  
  #randomregions_struc now contains:
  #the means of gene expresion (normalized with log, cutoff included)
  # per structure, per random region generated
  # so 1000 permutations = 1000 nested tables

  return(selected)
} 