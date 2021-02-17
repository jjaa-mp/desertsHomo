clean_perm_results <- function(permrun){
  #remove empty dataframes
  permrun <- permrun %>% 
    purrr::discard(empty)
  permrun <- lapply(permrun, as.data.frame)
  m_permrun <- melt(permrun)
  # m_permrun$variable <- stringr::str_remove_all(m_permrun$variable , "HSB.*[.]")
  m_permrun$L1 <- NULL
  m_permrun$variable <- as.character(m_permrun$variable)
  m_permrun$value <- as.numeric(m_permrun$value)
  m_permrun <- m_permrun %>% 
    dplyr::filter(value > 2)
  return(m_permrun)
} 
