#Internal function for permutation_test.R
get_data_perms <- function(mask) { 

  # Permutations
  randReg50set<-createRandomRegions(nregions = 1,length.mean = 15000000,
                                     length.sd = 1000000,genome = "hg19", mask = mask)
  #Changed region length 
  
  randomdesert <- c(runValue(seqnames(randReg50set)),start(randReg50set),end(randReg50set))
  # See comments in permutations code as to why this is implemented this way 
  chrstart <- paste0(randomdesert[1:2], collapse=":")
  result <- paste0(chrstart, "-",  randomdesert[3])
  return(result)

}