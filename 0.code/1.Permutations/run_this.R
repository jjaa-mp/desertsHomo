#load packages and functions for permutations
source("./packages.R")
lapply(list.files("./functions", full.names = TRUE), source)


# make sure to check file dependency paths- otherwise this will throw an error 
# load expression data - otherw
sestan = clean_sestan()
# load deserts of introgression data
akey = retrieve_GenesDeserts()


#Runs permutations in batches of 50
num <- 1:30
cleanperm <- NULL # contains no developmental stage data
rawperm <- NULL # contains developmental stage
for (i in num){
  permutationrun <- permutation_test(50, sestan, i)
  rawperm[[i]] <- permutationrun
  cleanperm[[i]] <- clean_perm_results(permutationrun)
}

# DESERTS
# reduces them to single permutation file
cleanperm <- cleanperm %>%  
  purrr::reduce(full_join)

#finally, run stat, plot differences, detect outliers
stats_sestperm_chen = stats_permutations(cleanperm, akey, sestan, "akey")
diffplot1 = difference_perm_sestan(stats_sestperm_chen, "sestan_", "chen", 8)
outliers_sestan_chen = detect_outlier_tb(stats_sestperm_chen, diffplot1, "sestan_chen")

#posthocs stages
rawperm <- lapply(rawperm, clean_rawperm)

rawperm <- rawperm %>% 
  purrr::reduce(full_join)

ph_chen_fetal1 = stats_permutations_stages(rawperm, akey, sestan, "akey", "fetal1")
ph_chen_fetal2 = stats_permutations_stages(rawperm, akey, sestan, "akey", "fetal2")
ph_chen_fetal3 = stats_permutations_stages(rawperm, akey, sestan, "akey", "fetal3")
ph_chen_birth_inf = stats_permutations_stages(rawperm, akey, sestan, "akey", "birth_inf")
ph_chen_inf_child = stats_permutations_stages(rawperm, akey, sestan, "akey", "inf_child")
ph_chen_child = stats_permutations_stages(rawperm, akey, sestan, "akey", "child")
ph_chen_adolescence = stats_permutations_stages(rawperm, akey, sestan, "akey", "adolescence")
ph_chen_adult = stats_permutations_stages(rawperm, akey, sestan, "akey", "adult")


# DESETS + POSITIVE SELECTION


stats_sestperm_chenpey = stats_permutations_s(permrun_s, pey_coords, sestan, "peycoords")
diffplot2 = difference_perm_sestan(stats_sestperm_chenpey, "sestan_", "chenpey", 8)
outliers_sestan_chenpey = detect_outlier_tb(stats_sestperm_chenpey, diffplot2, "sestan_chenpey")

