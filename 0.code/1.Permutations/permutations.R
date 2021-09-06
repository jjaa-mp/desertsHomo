
#Data dependencies: 
#- abadult
#- Sestan data
# Akey = Chen et al desert coordnitaes
# Pey = Peyrégne+deserts from Chen et al overlap

library(conflicted)
library(dotenv)
library(drake)
library(biomaRt)
library(ABAData)
library(ABAEnrichment)
library(plyr)
library(tidyverse)
library(GGally)
library(viridis)
library(ggpubr)
library(grid)
library(gridExtra)
library(lattice)
#library(xlsx)
library(readxl)
library(tidyr)
library(DescTools)
library(reshape2)
library(Matrix)
library(nlme)
library(multcomp)
library(ARTool)
library(emmeans)
library(PMCMRplus)
library(regioneR)
library(BSgenome.Hsapiens.UCSC.hg19)
library(R.utils)
library(rstatix)
library(outliers)
library(calibrate)

conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("reorder.factor", "gplots")
conflicted::conflict_prefer("melt", "reshape2")
conflicted::conflict_prefer("arrange", "dplyr")
conflicted::conflict_prefer("slice", "dplyr")
conflicted::conflict_prefer("strsplit", "base")
conflicted::conflict_prefer("select", "dplyr")

# --
# Functions:
# --

# Runs permutations:
permutation_test_s <- function(npermutations, dataset, seed) {
  #permutation_test, originally designed for Sestan data
  set.seed(seed) 
  
  mask <- data.frame(c("chr1","chr3","chr7","chr8"), 
                     c(105400000, 74100000,106200000,49400000), 
                     c(120600000, 89300000,123200000,66500000))
  
  #More efficient this way
  random <- replicate(n = npermutations, get_data_perms(mask), simplify = FALSE)
  run <- lapply(random, randomregion_biomart_query)
  #remove empty dataframes
  cleanrun <- run %>% 
    purrr::discard(empty)
  
  
  selected <- lapply(cleanrun, filter_genexpr_random, data = dataset)
  #filter_genexpr_random is an internal functions
  
  return(selected)
} 

#Internal function for permutation_test.R
filter_genexpr_random <- function(region, data) {
  subset.data <- data %>%
    dplyr::filter(gene_name %in% region$hgnc_symbol) %>% 
    dplyr::select(-c(gene_name))
  return(subset.data) 
}

# Plot 1
difference_perm_sestan <- function(stats_perm_df, where, which, height){
  
  stats_perm_df <- stats_perm_df %>% 
    group_by(struct_id, datasource) %>% 
    dplyr::summarize(mean_struct = mean(mean_struct))
  
  stats_perm_df$struct_id <- as.factor(stats_perm_df$struct_id)
  
  #stacked
  ggplot(stats_perm_df, aes (x= mean_struct, y = struct_id, fill = datasource )) +
    theme_minimal() + 
    geom_histogram(stat = "identity")
  
  diff <- stats_perm_df %>% 
    group_by(struct_id) %>% 
    dplyr::mutate(rest = mean_struct[datasource != 'permutations'] - mean_struct[datasource == 'permutations'] ) %>% 
    filter(datasource == 'permutations')
  
  #reorder
  diff$struct_id <- factor(diff$struct_id,                                    # Factor levels in decreasing order
                           levels = diff$struct_id[order(diff$rest, decreasing = FALSE)])
  
  a <-  ggplot(diff, aes (x= rest, y = struct_id )) +
    theme_minimal() + 
    labs(x = c("Log2 mean difference"), y = c("Structure")) +
    geom_histogram(stat = "identity",  color="steelblue", fill="steelblue") 
  
  ggsave(paste0("difference_perm_", where, which, ".pdf"), width = 8, height = height)
  
  return(diff)
}

# Statistics permutations: two similar functions per dataset, since there are diverging bits in the data structures
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
  
  return(full_data)
}
stats_permutations <- function(permutationrun, regionofinterest, aba, which){
  #Changes a bit compared to stats_permutations_s because of the ABA data structure
  
  perm_statsdf <- NULL
  perm_statsdf <- permutationrun
  perm_statsdf$datasource <- "permutations"
  colnames(perm_statsdf) <- c("idBrainStruct", "meansSubstruct", "datasource")
  
  #needed for later join
  perm_statsdf$idBrainStruct <- as.factor(perm_statsdf$idBrainStruct) 
  perm_statsdf$meansSubstruct <- as.numeric(perm_statsdf$meansSubstruct) 
  
  #checking normality
  #hist(perm_statsdf$meansSubstruct)
  #very normal looking, as expected
  
  #Gets mean expression from akey + sestan combination 
  aba_region_interest <- aba[aba$gene_name %in% regionofinterest$hgnc_symbol,]
  #aba_region_interest <- aba_region_interest[-1]
  #data wrangling for tidyness
  abakey_df <- melt(aba_region_interest)
  abakey_df <- abakey_df[-1]
  
  abakey_df <- abakey_df %>% 
    group_by(variable) %>% 
    dplyr::select(variable, value)  %>%  
    dplyr::filter(value > 2) %>% 
    dplyr::mutate(mean_struct = log2(value)) %>% 
    dplyr::select(variable, mean_struct) %>% 
    dplyr::summarize(mean_struct = mean(mean_struct)) %>% 
    dplyr::mutate(datasource = ifelse(which == "akey", "Chen", "Chen + Pey"))
  
  
  colnames(abakey_df) <- c("struct_id", "mean_struct", "datasource")
  
  #Checking normality
  colnames(perm_statsdf) <- c("struct_id", "mean_struct", "datasource")
  
  
  full_data <- full_join(perm_statsdf, abakey_df)
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
    write.csv(one.way2, paste0("output/ABAperm_posthoc_", "chen.csv"))
  } else if (which == "peycoords"){
    write.csv(one.way2, paste0("output/ABAperm_posthoc_", "chenpey.csv"))
  } 

  return(full_data)
}

# Posthocs per stages: sunukar to stats_permutations, but stage specific
stats_permutations_stages <- function(permrun_s, regionofinterest, sestan, which, stage){
  
  if(stage == "fetal1"){
    stag <- permrun_s[grep("153|150|113|103|149|114", permrun_s$variable), "variable"]
  } else if (stage == "fetal2") {
    stag <- permrun_s[grep("178|154|B96|B97", permrun_s$variable), "variable"]   
  } else if (stage == "fetal3") {
    stag <- permrun_s[grep("B98|107|B92|159", permrun_s$variable), "variable"]
  } else if (stage == "birth_inf") {
    stag <- permrun_s[grep("155|194|121|132|139", permrun_s$variable), "variable"]
  } else if (stage == "inf_child") {
    stag <- permrun_s[grep("131|171|122|143|173", permrun_s$variable), "variable"]
  } else if (stage == "child") {
    stag <- permrun_s[grep("172|118|141|174|175", permrun_s$variable), "variable"]
  } else if (stage == "adolescence") {
    stag <- permrun_s[grep("124|119|105|127", permrun_s$variable), "variable"]
  } else if (stage == "adult") {
    stag <- permrun_s[grep("130|136|126|145|123|135", permrun_s$variable), "variable"]
  }
  
  permrun_s <- permrun_s %>% 
    filter(variable %in% stag$variable)
  
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

# Outlier test
#diffdf = difference_perm_sestan returned object
detect_outlier_tb <- function(stats_dataframe, diffdf, which) {
  test <- stats_dataframe %>% 
    group_by(struct_id, datasource) %>% 
    dplyr::summarize(mean_struct = mean(mean_struct))
  
  test <- dcast(test,  struct_id ~ datasource)
  pdf(paste0(which, "_perm.pdf"))
  plot(test[[2]], test$permutations,
       pch = 16,
       xlab=which,
       ylab="permutations",
       textxy(test[[2]], test$permutations, test$struct_id, offset = -0.4))
  dev.off() #Saving the file
  

  
  a <- grubbs.test(diffdf$rest, opposite = FALSE)
  b <- grubbs.test(diffdf$rest, opposite = TRUE)
  
  Matrix::print(diffdf[which.max(diffdf$rest),])
  Matrix::print(a)
  
  Matrix::print(diffdf[which.min(diffdf$rest),])
  Matrix::print(b)
  
}


# --
# Permutations for ABA data
# --

permutationrunABA1 = permutation_test_s(50,abadult,100) # 100 = the seed; has to change from perm run to perm run
cleanpermABA1 = clean_perm_results_s(permutationrunABA1)
permutationrunABA2 = permutation_test_s(50,abadult,101)
cleanpermABA2 = clean_perm_results_s(permutationrunABA2)
permutationrunABA3 = permutation_test_s(50,abadult,103)
cleanpermABA3 = clean_perm_results_s(permutationrunABA3)
permutationrunABA4 = permutation_test_s(50,abadult,104)
cleanpermABA4 = clean_perm_results_s(permutationrunABA4)
permutationrunABA5 = permutation_test_s(50,abadult,105)
cleanpermABA5 = clean_perm_results_s(permutationrunABA5)
permutationrunABA6 = permutation_test_s(50,abadult,106)
cleanpermABA6 = clean_perm_results_s(permutationrunABA6)
permutationrunABA7 = permutation_test_s(50,abadult,107)
cleanpermABA7 = clean_perm_results_s(permutationrunABA7)
permutationrunABA8 = permutation_test_s(50,abadult,108)
cleanpermABA8 = clean_perm_results_s(permutationrunABA8)
permutationrunABA9 = permutation_test_s(50,abadult,109)
cleanpermABA9 = clean_perm_results_s(permutationrunABA9)
permutationrunABA10 = permutation_test_s(50,abadult,110)
cleanpermABA10 = clean_perm_results_s(permutationrunABA10) #cleans, normalizes by log2
permutationrun = rbind(cleanpermABA1, cleanpermABA2, #joins the results int
                       cleanpermABA3, cleanpermABA4,
                       cleanpermABA5, cleanpermABA6,
                       cleanpermABA7, cleanpermABA8,
                       cleanpermABA9, cleanpermABA10)

#Statistics, including posthoc, + difference plot  outlier tests
permutationstats = stats_permutations(permutationrun, akey, abadult, "akey")
pdiffplot1 = difference_perm_sestan(permutationstats, "ABA_",  "chen", 40)
outliers_ABA_chen = detect_outlier_tb(permutationstats, pdiffplot1, "ABA_chen")

permutationstats_pey = stats_permutations(permutationrun, pey_coords, abadult, "peycoords")
pdiffplot2 = difference_perm_sestan(permutationstats_pey, "ABA_", "chenpey", 40)
outliers_ABA_chenpey = detect_outlier_tb(permutationstats_pey, pdiffplot2, "ABA_chenpey")

# --
# Permutations for Sestan data
# --

sestan = clean_sestan()
permutationrun1 = permutation_test_s(50, sestan, 100) # 100 = the seed has to change from perm run to perm run
cleanperm1 = clean_perm_results(permutationrun1)
permutationrun2 = permutation_test_s(50, sestan, 101)
cleanperm2 = clean_perm_results(permutationrun2)
permutationrun3 = permutation_test_s(50, sestan, 102)
cleanperm3 = clean_perm_results(permutationrun3)
permutationrun4 = permutation_test_s(50, sestan, 103)
cleanperm4 = clean_perm_results(permutationrun4)
permutationrun5 = permutation_test_s(50, sestan, 104)
cleanperm5 = clean_perm_results(permutationrun5)
permutationrun6 = permutation_test_s(50, sestan, 105)
cleanperm6 = clean_perm_results(permutationrun6)
permutationrun7 = permutation_test_s(50, sestan, 106)
cleanperm7 = clean_perm_results(permutationrun7)
permutationrun8 = permutation_test_s(50, sestan, 107)
cleanperm8 = clean_perm_results(permutationrun8)
permutationrun9 = permutation_test_s(50, sestan, 108)
cleanperm9 = clean_perm_results(permutationrun9)
permutationrun10 = permutation_test_s(50, sestan, 109)
cleanperm10 = clean_perm_results(permutationrun10)  #cleans, normalizes by log2
permrun_s = rbind(cleanperm1, cleanperm2, 
                  cleanperm3, cleanperm4, 
                  cleanperm5, cleanperm6, 
                  cleanperm7, cleanperm8, 
                  cleanperm9, cleanperm10)

#Statistics, including posthoc, + difference plot  outlier tests
stats_sestperm_chen = stats_permutations(permrun_s, akey, sestan, "akey")
diffplot1 = difference_perm_sestan(stats_sestperm_chen, "sestan_", "chen", 8)
outliers_sestan_chen = detect_outlier_tb(stats_sestperm_chen, diffplot1, "sestan_chen")

stats_sestperm_chenpey = stats_permutations_s(permrun_s, pey_coords, sestan, "peycoords")
diffplot2 = difference_perm_sestan(stats_sestperm_chenpey, "sestan_", "chenpey", 8)
outliers_sestan_chenpey = detect_outlier_tb(stats_sestperm_chenpey, diffplot2, "sestan_chenpey")

#Posthocs for specific stages
ph_chen_fetal1 = stats_permutations_stages(rawperm, akey, sestan, "akey", "fetal1")
ph_chen_fetal2 = stats_permutations_stages(permrun_s, akey, sestan, "akey", "fetal2")
ph_chen_fetal3 = stats_permutations_stages(permrun_s, akey, sestan, "akey", "fetal3")
ph_chen_birth_inf = stats_permutations_stages(permrun_s, akey, sestan, "akey", "birth_inf")
ph_chen_inf_child = stats_permutations_stages(permrun_s, akey, sestan, "akey", "inf_child")
ph_chen_child = stats_permutations_stages(permrun_s, akey, sestan, "akey", "child")
ph_chen_adolescence = stats_permutations_stages(permrun_s, akey, sestan, "akey", "adolescence")
ph_chen_adult = stats_permutations_stages(permrun_s, akey, sestan, "akey", "adult")

ph_chenpey_fetal1 = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "fetal1")
ph_chenpey_fetal2 = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "fetal2")
ph_chenpey_fetal3 = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "fetal3")
ph_chenpey_birth_inf = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "birth_inf")
ph_chenpey_inf_child = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "inf_child")
ph_chenpey_child = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "child")
ph_chenpey_adolescence = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "adolescence")
ph_chenpey_adult = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "adult")

