stats_permutations_s <- function(permrun_s, regionofinterest, sestan, which){
  perm_statsdf <- NULL
  perm_statsdf <- permrun_s
  perm_statsdf$datasource <- "permutations"
  colnames(perm_statsdf) <- c("idBrainStruct", "meansSubstruct", "datasource")
  
  #needed for later join
  perm_statsdf$idBrainStruct <- as.factor(perm_statsdf$idBrainStruct) 
  perm_statsdf$meansSubstruct <- as.numeric(perm_statsdf$meansSubstruct) 
  
  #checking normality
  #hist(perm_statsdf$meansSubstruct)
  #very normal looking, as expected
  
  #Gets mean expression from akey + sestan combination 
  sestan_region_interest <- sestan[sestan$gene_name %in% regionofinterest$hgnc_symbol,]
  sestan_region_interest <- sestan_region_interest[-1]
  #data wrangling for tidyness
  sestakey_df <- melt(sestan_region_interest)
  
  sestakey_df <- sestakey_df %>% 
    group_by(variable) %>% 
    dplyr::select(variable, value)  %>%  
    dplyr::filter(value > 2) %>% 
    dplyr::mutate(mean_struct = log2(value)) %>% 
    dplyr::select(variable, mean_struct) %>% 
    dplyr::summarize(mean_struct = mean(mean_struct)) %>% 
    dplyr::mutate(datasource = ifelse(which == "akey", "Chen", "Chen + Pey"))
  
  
  colnames(sestakey_df) <- c("struct_id", "mean_struct", "datasource")
  
  #Checking normality
  #hist(abakey_df$mean_struct) 
  # Enough
  colnames(perm_statsdf) <- c("struct_id", "mean_struct", "datasource")
  
  
  full_data <- full_join(perm_statsdf, sestakey_df)
  full_data$datasource <- as.factor(full_data$datasource)
  #quick visualization
  ggplot(full_data) +
    theme_minimal() +
    aes(x = mean_struct, color = datasource, fill = datasource) +
    geom_histogram() +
    theme(legend.position = "top")
  
  full_data$struct_id <- stringr::str_remove_all(full_data$struct_id, ".*[.]")
  
  fulldata <- full_data %>% 
    group_by(datasource, struct_id) %>% 
    dplyr::summarize(mean_struct = mean(mean_struct)) # since duplicated ids happen 
  #this differs from regular stats_permputations script
  
  #from car package: tests if variances are equal (they are)
  #leveneTest(mean_struct ~ struct_id, data = full_data)
  
  printf("Wait a bit; at 1000 permutations, this should take some time")
  res.aov <- anova_test(data = full_data, dv = mean_struct, wid = struct_id, between = datasource)
  res.aov
  printf("Done!")
  
  full_data$struct_id <- as.character(full_data$struct_id)
  printf("Posthoc by structure")
  one.way2 <- full_data %>%
    group_by(struct_id) %>%
    anova_test(dv = mean_struct, between  = datasource) %>%
    get_anova_table() %>%
    adjust_pvalue(method = "bonferroni")
  
  one.way2

  if (which == "akey"){
    write.csv(one.way2, paste0("output/sestanperm_posthoc_", "chen.csv"))
  } else if (which == "peycoords"){
    write.csv(one.way2, paste0("output/sestanperm_posthoc_", "chenpey.csv"))
  } 
  
  
  #if (which == "akey"){
  #  plotting_topbottom(full_data, "sestantop", "chen")
  #  plotting_topbottom(full_data, "sestanbottom", "chen")
  #} else if (which == "peycoords"){
  #  plotting_topbottom(full_data, "sestantop", "chenpey")
  #  plotting_topbottom(full_data, "sestanbottom", "chenpey")
  #} 
  return(full_data)
}