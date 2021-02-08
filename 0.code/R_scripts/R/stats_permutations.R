stats_permutations <- function(permutationrun, regionofinterest, abdata, which){
  means <- lapply(permutationrun, as.data.frame)
  names <- means[[1]]$idBrainStruct
  means <- lapply(means, '[[', 2)
  means <- unlist(means)
  #All the means of 414 organs
  
  perm_statsdf <- NULL
  perm_statsdf$struct_id <- names
  perm_statsdf$mean_struct <- means
  perm_statsdf <- as.data.frame(perm_statsdf)
  perm_statsdf$datasource <- "permutations"
  
  #needed for later join
  perm_statsdf$struct_id <- as.factor(perm_statsdf$struct_id) 
  perm_statsdf$mean_struct <- as.numeric(perm_statsdf$mean_struct) 
  
  #checking normality
  #hist(perm_statsdf$mean_struct)
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
    dplyr::select(variable, mean_struct) %>% 
    dplyr::mutate(datasource = ifelse(which == "akey", "Chen", "Chen + Pey")) 
  
  colnames(abakey_df) <- c("struct_id", "mean_struct", "datasource")

  #Checking normality
  #hist(abakey_df$mean_struct) 
  # Enough
  
  full_data <- full_join(perm_statsdf, abakey_df)
  
  #quick visualization
  ggplot(full_data) +
    theme_minimal() +
    aes(x = mean_struct, color = datasource, fill = datasource) +
    geom_histogram() +
    theme(legend.position = "none")
  
  #from car package: tests if variances are equal (they are)
  #leveneTest(mean_struct ~ struct_id, data = full_data)
  
  printf("Wait a bit; at 1000 permutations, this should take some time")
  res.aov <- anova_test(data = full_data, dv = mean_struct, wid = struct_id, between = datasource)
  #printf(res.aov)
  printf("Done!")
  
  printf("Posthoc by structure")
  one.way2 <- full_data %>%
    group_by(struct_id) %>%
    anova_test(dv = mean_struct, between  = datasource) %>%
    get_anova_table() %>%
    adjust_pvalue(method = "bonferroni")
  
  one.way2
  #printf(one.way2)
  if (which == "akey"){
    write.csv(one.way2, paste0("output/perm_posthoc_", "chen.csv"))
  } else if (which == "peycoords"){
    write.csv(one.way2, paste0("output/perm_posthoc_", "chenpey.csv"))
  } 

  
  if (which == "akey"){
    plotting_topbottom(full_data, "top", "chen")
    plotting_topbottom(full_data, "bottom", "chen")
  } else if (which == "peycoords"){
    plotting_topbottom(full_data, "top", "chenpey")
    plotting_topbottom(full_data, "bottom", "chenpey")
  } 
  return(full_data)

  #residuals, with some artifacts from log transformation I guess
  #plot(fitted(model),
  #           residuals(model))
}



