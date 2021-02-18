# list of functions to call, as stored in .R files

the_plan <- drake_plan(
  	akey = retrieve_GenesDeserts(),
  	abadult = dt_adults(), #gets ab data
  	
  	pey = load_data("file_dependencies/2020_pey_coords.bed"), 
  	pey_coords = biominput(pey, akey, TRUE), # Genes within Akey and Peyregne
  	rac = load_data("file_dependencies/2020_rac_coords.bed"),
  	rac_coords = biominput(rac, akey, TRUE), # Genes within Akey and Racimo

  	abadult_pl = abadult_plots(abadult,  #gets ab plots
  	              akey, 
  	              pey_coords),
  	ab = abad_5_stages(), 
  	q75(ab, akey),  #gets quantile 75 plots in deserts
  	
  	#The following block should be checked to introduce:
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
  # 
  # 
  #  #Permutations  (using ABA data)
  #npermutations=1000,
  #permutationrun = permutation_test(npermutations, abadult),
  #permutationstats = stats_permutations(permutationrun, akey, abadult, "akey"),
  #permutationstats_pey = stats_permutations(permutationrun, pey_coords, abadult, "peycoords"),
    
  #sestan_a = mRNA_sestan(akey),
  #sestan_apey = mRNA_sestan(pey_coords),
  #sestan_raw = mRNA_sestan() #Done in cluster
  #combinedplot_sestan = plot_sestan_compared(sestan_a, sestan_apey),
  
  #Permutations (using Sestan data)
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
  cleanperm10 = clean_perm_results_s(permutationrun10),
  permrun_s = rbind(cleanperm1, 
                    cleanperm2, 
                    cleanperm3, 
                    cleanperm4, 
                    cleanperm5, 
                    cleanperm6, 
                    cleanperm7, 
                    cleanperm8, 
                    cleanperm9, 
                    cleanperm10),
  
  stats_sestperm_chen = stats_permutations_s(permrun_s, akey, sestan, "akey"),
  diffplot1 = difference_perm_sestan(stats_sestperm_chen, "chen"),
  stats_sestperm_chenpey = stats_permutations_s(permrun_s, pey_coords, sestan, "peycoords"),
  diffplot2 = difference_perm_sestan(stats_sestperm_chenpey, "chenpey")
)

