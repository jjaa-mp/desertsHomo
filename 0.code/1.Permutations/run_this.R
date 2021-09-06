#load packages and functions for permutations
source("./packages.R")
lapply(list.files("./functions", full.names = TRUE), source)


# --
# Dependencies
# --

# NOTE:
# make sure to check file path (ideally, provide them in a folder in the
# working directiory called "file dependencies/" with the, well, dependencies
# - otherwise this will throw an error 

# load expression data for sestan and the adult stages of ABA
sestan = clean_sestan()
abadult = dt_adults()
# load deserts of introgression data
akey = retrieve_GenesDeserts()


# --
# Permutations for Sestan data
# --


#Runs permutations in batches of 50
num <- 1:30
cleanperm <- NULL # contains no developmental stage data
rawperm <- NULL # contains developmental stage
for (i in num){
  permutationrun <- permutation_test(50, sestan, i)
  rawperm[[i]] <- permutationrun
  cleanperm[[i]] <- clean_perm_results(permutationrun)
}

# SESTAN: DESERTS
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


# SESTAN: DESETS + POSITIVE SELECTION
pey <- read_table2("file_dependencies/2020_pey_coords.bed", 
                                col_names = FALSE)
pey_coords <- biominput(pey, akey, TRUE) # Genes within Akey and Peyregne


stats_sestperm_chenpey = stats_permutations(cleanperm, pey_coords, sestan, "peycoords")
diffplot2 = difference_perm_sestan(stats_sestperm_chenpey, "sestan_", "chenpey", 8)
outliers_sestan_chenpey = detect_outlier_tb(stats_sestperm_chenpey, diffplot2, "sestan_chenpey")


# posthocs per stage
ph_chenpey_fetal1 = stats_permutations_stages(rawperm, pey_coords, sestan, "peycoords", "fetal1")
ph_chenpey_fetal2 = stats_permutations_stages(rawperm, pey_coords, sestan, "peycoords", "fetal2")
ph_chenpey_fetal3 = stats_permutations_stages(rawperm, pey_coords, sestan, "peycoords", "fetal3")
ph_chenpey_birth_inf = stats_permutations_stages(rawperm, pey_coords, sestan, "peycoords", "birth_inf")
ph_chenpey_inf_child = stats_permutations_stages(rawperm, pey_coords, sestan, "peycoords", "inf_child")
ph_chenpey_child = stats_permutations_stages(rawperm, pey_coords, sestan, "peycoords", "child")
ph_chenpey_adolescence = stats_permutations_stages(rawperm, pey_coords, sestan, "peycoords", "adolescence")
ph_chenpey_adult = stats_permutations_stages(rawperm, pey_coords, sestan, "peycoords", "adult")



# --
# Permutations for ABA data
# --

#Runs permutations in batches of 50
num <- 1:30
cleanperm <- NULL 
rawperm <- NULL 
for (i in num){
  permutationrun <- permutation_test(50, abadult, i)
  rawperm[[i]] <- permutationrun
  cleanperm[[i]] <- clean_perm_results(permutationrun)
}

# reduces them to single permutation file
cleanperm <- cleanperm %>%  
  purrr::reduce(full_join)

# ABA: DESERTS
permutationstats = stats_permutations(cleanperm, akey, abadult, "akey")
pdiffplot1 = difference_perm_sestan(permutationstats, "ABA_",  "chen", 40)
outliers_ABA_chen = detect_outlier_tb(permutationstats, pdiffplot1, "ABA_chen")

permutationstats_pey = stats_permutations(cleanperm, pey_coords, abadult, "peycoords")
pdiffplot2 = difference_perm_sestan(permutationstats_pey, "ABA_", "chenpey", 40)
outliers_ABA_chenpey = detect_outlier_tb(permutationstats_pey, pdiffplot2, "ABA_chenpey")
