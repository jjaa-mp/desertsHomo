clean_perm_results_s <- function(permrun){
  #clean, but for sestan
  permrun <- melt(permrun)
  permrun <- permrun %>% 
    group_by(variable) %>%
    filter(value > 2)
  permrun$value <- log2(permrun$value) 
  permrun <- permrun %>% 
    dplyr::summarize(value = mean(value))
}


