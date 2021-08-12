stats_permutations_stages <- function(permrun_s, regionofinterest, sestan, which, stage){
  
  if(stage == "fetal1"){
    stag$idBrainStruct <- permrun_s[grep("153|150|113|103|149|114", permrun_s$idBrainStruct), "idBrainStruct"]
  } else if (stage == "fetal2") {
    stag$idBrainStruct <- permrun_s[grep("178|154|B96|B97", permrun_s$idBrainStruct), "idBrainStruct"]   
  } else if (stage == "fetal3") {
    stag$idBrainStruct <- permrun_s[grep("B98|107|B92|159", permrun_s$idBrainStruct), "idBrainStruct"]
  } else if (stage == "birth_inf") {
    stag$idBrainStruct <- permrun_s[grep("155|194|121|132|139", permrun_s$idBrainStruct), "idBrainStruct"]
  } else if (stage == "inf_child") {
    stag$idBrainStruct <- permrun_s[grep("131|171|122|143|173", permrun_s$idBrainStruct), "idBrainStruct"]
  } else if (stage == "child") {
    stag$idBrainStruct <- permrun_s[grep("172|118|141|174|175", permrun_s$idBrainStruct), "idBrainStruct"]
  } else if (stage == "adolescence") {
    stag$idBrainStruct <- permrun_s[grep("124|119|105|127", permrun_s$idBrainStruct), "idBrainStruct"]
  } else if (stage == "adult") {
    stag$idBrainStruct <- permrun_s[grep("130|136|126|145|123|135", permrun_s$idBrainStruct), "idBrainStruct"]
  }
  
  permrun_s <- permrun_s %>% 
    filter(idBrainStruct %in% stag$idBrainStruct)
  
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
  
  if(stage == "fetal1"){
    stag <- sestakey_df[grep("153|150|113|103|149|114", sestakey_df$struct_id), "struct_id"]
  } else if (stage == "fetal2") {
    stag <- sestakey_df[grep("178|154|B96|B97", sestakey_df$struct_id), "struct_id"]   
  } else if (stage == "fetal3") {
    stag <- sestakey_df[grep("B98|107|B92|159", sestakey_df$struct_id), "struct_id"]
  } else if (stage == "birth_inf") {
    stag <- sestakey_df[grep("155|194|121|132|139", sestakey_df$struct_id), "struct_id"]
  } else if (stage == "inf_child") {
    stag <- sestakey_df[grep("131|171|122|143|173", sestakey_df$struct_id), "struct_id"]
  } else if (stage == "child") {
    stag <- sestakey_df[grep("172|118|141|174|175", sestakey_df$struct_id), "struct_id"]
  } else if (stage == "adolescence") {
    stag <- sestakey_df[grep("124|119|105|127", sestakey_df$struct_id), "struct_id"]
  } else if (stage == "adult") {
    stag <- sestakey_df[grep("130|136|126|145|123|135", sestakey_df$struct_id), "struct_id"]
  }
  
  
  sestakey_df <- sestakey_df %>% 
    filter(struct_id %in% stag$struct_id)
  
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
    write.csv(one.way2, paste0("output/sestanperm_ph_", stage, "_chen.csv"))
  } else if (which == "peycoords"){
    write.csv(one.way2, paste0("output/sestanperm_ph_", stage, "_chenpey.csv"))
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
