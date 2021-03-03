# list of functions to call, as stored in .R files

the_plan <- drake_plan(
  	akey = retrieve_GenesDeserts(),
  	abadult = dt_adults(), #gets ab data
  	
  	pey = load_data("file_dependencies/2020_pey_coords.bed"), 
  	pey_coords = biominput(pey, akey, TRUE), # Genes within Akey and Peyregne
  	rac = load_data("file_dependencies/2020_rac_coords.bed"), #  Genes within Racimo
  	rac_coords = biominput(rac, akey, TRUE), # Genes within Akey and Racimo

  	abadult_pl = abadult_plots(abadult,  #gets ab plots
  	              akey, 
  	              pey_coords),
  	ab = abad_5_stages(), 
  	q75(ab, akey),  #gets quantile 75 plots in deserts
  
  #  #Permutations  (using ABA data)
  #Currently has to run in 10 separate runs
  permutationrunABA1 = permutation_test_s(100,abadult,100), # 100 = the seed; has to change from perm run to perm run
  cleanpermABA1 = clean_perm_results_s(permutationrunABA1),
  permutationrunABA2 = permutation_test_s(100,abadult,101),
  cleanpermABA2 = clean_perm_results_s(permutationrunABA2),
  permutationrunABA3 = permutation_test_s(100,abadult,103),
  cleanpermABA3 = clean_perm_results_s(permutationrunABA3),
  permutationrunABA4 = permutation_test_s(100,abadult,104),
  cleanpermABA4 = clean_perm_results_s(permutationrunABA4),
  permutationrunABA5 = permutation_test_s(100,abadult,105),
  cleanpermABA5 = clean_perm_results_s(permutationrunABA5),
  permutationrunABA6 = permutation_test_s(100,abadult,106),
  cleanpermABA6 = clean_perm_results_s(permutationrunABA6),
  permutationrunABA7 = permutation_test_s(100,abadult,107),
  cleanpermABA7 = clean_perm_results_s(permutationrunABA7),
  permutationrunABA8 = permutation_test_s(100,abadult,108),
  cleanpermABA8 = clean_perm_results_s(permutationrunABA8),
  permutationrunABA9 = permutation_test_s(100,abadult,100),
  cleanpermABA9 = clean_perm_results_s(permutationrunABA9),
  permutationrunABA10 = permutation_test_s(100,abadult,100),
  cleanpermABA10 = clean_perm_results_s(permutationrunABA10), #cleans, normalizes by log2
  permutationrun = rbind(cleanpermABA1, cleanpermABA2, #joins the results int
                         cleanpermABA3, cleanpermABA4,
                         cleanpermABA5, cleanpermABA6,
                         cleanpermABA7, cleanpermABA8,
                         cleanpermABA9, cleanpermABA10),
  
  #Statistics, including posthoc, + difference plot  outlier tests
  permutationstats = stats_permutations(permutationrun, akey, abadult, "akey"),
  pdiffplot1 = difference_perm_sestan(permutationstats, "ABA_",  "chen", 40),
  outliers_ABA_chen = detect_outlier_tb(permutationstats, pdiffplot1, "ABA_chen"),
  
  permutationstats_pey = stats_permutations(permutationrun, pey_coords, abadult, "peycoords"),
  pdiffplot2 = difference_perm_sestan(permutationstats_pey, "ABA_", "chenpey", 40),
  outliers_ABA_chenpey = detect_outlier_tb(permutationstats_pey, pdiffplot2, "ABA_chenpey"),
  

  #Permutations (using Sestan data)
  #Currently has to run in 10 separate runs
  sestan = clean_sestan(),
  permutationrun1 = permutation_test_s(100, sestan, 100), # 100 = the seed has to change from perm run to perm run
  cleanperm1 = clean_perm_results_s(permutationrun1),
  permutationrun2 = permutation_test_s(100, sestan, 101),
  cleanperm2 = clean_perm_results_s(permutationrun2),
  permutationrun3 = permutation_test_s(100, sestan, 102),
  cleanperm3 = clean_perm_results_s(permutationrun3),
  permutationrun4 = permutation_test_s(100, sestan, 103),
  cleanperm4 = clean_perm_results_s(permutationrun4),
  permutationrun5 = permutation_test_s(100, sestan, 104),
  cleanperm5 = clean_perm_results_s(permutationrun5),
  permutationrun6 = permutation_test_s(100, sestan, 105),
  cleanperm6 = clean_perm_results_s(permutationrun6),
  permutationrun7 = permutation_test_s(100, sestan, 106),
  cleanperm7 = clean_perm_results_s(permutationrun7),
  permutationrun8 = permutation_test_s(100, sestan, 107),
  cleanperm8 = clean_perm_results_s(permutationrun8),
  permutationrun9 = permutation_test_s(100, sestan, 108),
  cleanperm9 = clean_perm_results_s(permutationrun9),
  permutationrun10 = permutation_test_s(100, sestan, 109),
  cleanperm10 = clean_perm_results_s(permutationrun10),  #cleans, normalizes by log2
  permrun_s = rbind(cleanperm1, cleanperm2, 
                    cleanperm3, cleanperm4, 
                    cleanperm5, cleanperm6, 
                    cleanperm7, cleanperm8, 
                    cleanperm9, cleanperm10),
  
  #Statistics, including posthoc, + difference plot  outlier tests
  stats_sestperm_chen = stats_permutations_s(permrun_s, akey, sestan, "akey"),
  diffplot1 = difference_perm_sestan(stats_sestperm_chen, "sestan_", "chen", 8),
  outliers_sestan_chen = detect_outlier_tb(stats_sestperm_chen, diffplot1, "sestan_chen"),
  
  stats_sestperm_chenpey = stats_permutations_s(permrun_s, pey_coords, sestan, "peycoords"),
  diffplot2 = difference_perm_sestan(stats_sestperm_chenpey, "sestan_", "chenpey", 8),
  outliers_sestan_chenpey = detect_outlier_tb(stats_sestperm_chenpey, diffplot2, "sestan_chenpey"),

  #Posthocs for specific stages
  ph_chen_fetal1 = stats_permutations_stages(permrun_s, akey, sestan, "akey", "fetal1"),
  ph_chen_fetal2 = stats_permutations_stages(permrun_s, akey, sestan, "akey", "fetal2"),
  ph_chen_fetal3 = stats_permutations_stages(permrun_s, akey, sestan, "akey", "fetal3"),
  ph_chen_birth_inf = stats_permutations_stages(permrun_s, akey, sestan, "akey", "birth_inf"),
  ph_chen_inf_child = stats_permutations_stages(permrun_s, akey, sestan, "akey", "inf_child"),
  ph_chen_child = stats_permutations_stages(permrun_s, akey, sestan, "akey", "child"),
  ph_chen_adolescence = stats_permutations_stages(permrun_s, akey, sestan, "akey", "adolescence"),
  ph_chen_adult = stats_permutations_stages(permrun_s, akey, sestan, "akey", "adult"),
  
  ph_chenpey_fetal1 = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "fetal1"),
  ph_chenpey_fetal2 = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "fetal2"),
  ph_chenpey_fetal3 = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "fetal3"),
  ph_chenpey_birth_inf = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "birth_inf"),
  ph_chenpey_inf_child = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "inf_child"),
  ph_chenpey_child = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "child"),
  ph_chenpey_adolescence = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "adolescence"),
  ph_chenpey_adult = stats_permutations_stages(permrun_s, pey_coords, sestan, "peycoords", "adult"),

  
  orderperm1 = permrun_order(permutationrun1),
  orderperm2 =  permrun_order(permutationrun2),
  orderperm3 = permrun_order(permutationrun3),
  orderperm4 = permrun_order(permutationrun4),
  orderperm5 = permrun_order(permutationrun5),
  orderperm6 = permrun_order(permutationrun6),
  orderperm7 = permrun_order(permutationrun7),
  orderperm8 = permrun_order(permutationrun8),
  orderperm9 = permrun_order(permutationrun9),
  orderperm10 = permrun_order(permutationrun10),
  
  ordered_all = allperm_order(orderperm1, orderperm2,
                              orderperm3, orderperm4,
                              orderperm5, orderperm6,
                              orderperm7, orderperm8,
                              orderperm9, orderperm10),
  likelihood_chen = top_structures_how_likely(ordered_all, stats_sestperm_chen, "Deserts of introgression"),
  likelihood_chenpey = top_structures_how_likely(ordered_all, stats_sestperm_chenpey, "Positively selected in deserts of introgression"),
  
)


#The following blocks should be checked to introduce:
#-the new log normalization 
#-erase non-parametric tests

# 	abakey_data = ab_data_plots(ab, akey, TRUE), # prepares akey data for ab plots
# 	abakey_plots = ab_plots(abakey_data, "Akey", " Deserts", "Striatum"),  #generates plots
# 	
# 	abakeypey_data = ab_data_plots(ab, pey_coords, TRUE), # prepares pey+akey data for ab plots
# 	abakepey_plots = ab_plots(abakeypey_data, "AkeyPey", " Deserts and Pey", "Somato - Motor - Parietal - Aud Ctxm"),
# 	generates plots
# 	
# 	abakeyrac_data = ab_data_plots(ab, rac_coords, TRUE), # prepares rac+akey data for ab plots
# 	abakeyrac_plots = ab_plots(abakeyrac_data, "AkeyRac", " Rac and Akey", "Striatum"),  #generates plots
# 	
# 	rac_coords_noakey = biominput(rac, akey, FALSE), # prepares rac data for ab plots
# 	abrac_data = ab_data_plots(ab, rac_coords_noakey, TRUE),
# 	abrac_plots = ab_plots(abrac_data, "Rac", " Rac", "Striatum"), #generates plots
# 	
# 	pey_coords_noakey = biominput(pey, akey, FALSE), # prepares pey data for ab plots
# 	abpey_data = ab_data_plots(ab, pey_coords_noakey, TRUE),
# 	abpey_plots = ab_plots(abpey_data, "Pey", " Pey", "Striatum"),
# 	 
# 	abakeypey_stats = ab_data_plots(ab, pey_coords, FALSE), #without means to calculate stats
# 	anova_akeypey = test_nonparam_anova(abakeypey_stats, "chenpey"), 
# 	fr_akeypey = friedman_only(abakeypey_stats, "chenpey"),
# 	nemenyi_posthoc(abakeypey_stats, "chenpey"),
# 	
# 	abakey_stats = ab_data_plots(ab, akey, FALSE),  #without means to calculate stats
# 	anova_akey = test_nonparam_anova(abakey_stats, "chen"),
# 	fr_akey = friedman_only(abakey_stats, "chen"),
# 	nemenyi_posthoc(abakey_stats, "chen"),
# 	 
# 	GSEakey = GSE_data(akey, pey_coords, "akey"), # prepares GSE data
# 	ttest_result = ttest(GSEakey), 
# 	akey_fried = friedman(GSEakey, "akey"), # friedman test for akey
# 	#anova_GSEakey = alt_anova(GSEakey), # checking with anova + Tukey if Friedman is alright
# 	 
# 	GSEpey = GSE_data(akey, pey_coords, "akeypey"),
# 	pey_fried = friedman(GSEpey, "pey"), # friedman test for pey
#  #anova_GSEpey = alt_anova(GSEpey) # checking with anova + Tukey if Friedman is alright

