quantilecut <- function(df) {
  cutoff <- df %>% dplyr::filter(. > quantile(0.10))
  return(cutoff)
} 
