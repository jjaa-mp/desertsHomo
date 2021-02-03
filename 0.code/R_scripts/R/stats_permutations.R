stats_permutations <- function(permutationrun, regionofinterest, abdata){
  means <- lapply(permutationrun, as.data.frame)
  names <- means[[1]]$idBrainStruct
  means <- lapply(means, '[[', 2)
  means <- unlist(means)
  #All the means of 414 organs
  
  #tidy version of nested list (permutationrun)
  perm_statsdf <- NULL
  perm_statsdf$struct_id <- names
  perm_statsdf$mean_struct <- means
  perm_statsdf <- as.data.frame(perm_statsdf)
  
  #needed for later join
  perm_statsdf$struct_id <- as.factor(perm_statsdf$struct_id) 
  perm_statsdf$mean_struct <- as.numeric(perm_statsdf$mean_struct) 
  
  #checking normality
  hist(perm_statsdf$mean_struct)
  #very normal looking, as expected
  

  #Gets mean expression from akey + brain allen (adults) combination 
  ab_region_interest <- abdata[abdata$gene_name %in% regionofinterest$hgnc_symbol,]
  
  #data wrangling for tidyness
  abakey_df <- melt(ab_region_interest)
  
  abakey_df <- abakey_df %>% 
    group_by(variable) %>% 
    dplyr::select(variable, value)  %>% 
    dplyr::filter(value > quantile(value, 0.10)) %>% 
    dplyr::mutate(mean_struct = log(value)) %>% 
    dplyr::select(variable, mean_struct) 
  
  colnames(abakey_df) <- c("struct_id", "mean_struct")

  #Checking normality
  hist(abakey_df$mean_struct) 
  # Enough
  
  test <- full_join(perm_statsdf, abakey_df)
}




