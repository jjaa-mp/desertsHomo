clean_perm_results <- function(permrun){
  #remove empty dataframes
  permrun <- permrun %>% 
    purrr::discard(empty)
  permrun <- lapply(permrun, as.data.frame)
  m_permrun <- reshape2::melt(permrun)
  m_permrun$idBrainStruct <- stringr::str_remove_all(m_permrun$idBrainStruct , "HSB.*[.]")
  m_permrun$L1 <- NULL
  m_permrun$idBrainStruct <- as.character(m_permrun$idBrainStruct)
  m_permrun$meansSubstruct <- as.numeric(m_permrun$meansSubstruct)
  return(m_permrun)
} 
