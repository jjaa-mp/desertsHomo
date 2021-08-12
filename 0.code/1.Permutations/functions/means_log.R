#Internal to permutation_test.R
means_log <- function(abDeserts){
  meansSubstruct <- colMeans(log2(abDeserts))
  rownames(meansSubstruct)
  meansSubstructTemp <- cbind(idBrainStruct = rownames(meansSubstruct), meansSubstruct)
  meansSubstructTemp <- cbind(idBrainStruct = rownames(meansSubstructTemp), meansSubstructTemp)
  rownames(meansSubstructTemp) <- 1:nrow(meansSubstructTemp)
  return(meansSubstructTemp)
}