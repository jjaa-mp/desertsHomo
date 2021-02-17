permutation_test_s <- function(npermutations, dataset, seed) {
  #permutation_test, but for Sestan data
  set.seed(seed) # won't work if lock_envir = TRUE in _drake.R
  
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

  
  return(selected)
} 
