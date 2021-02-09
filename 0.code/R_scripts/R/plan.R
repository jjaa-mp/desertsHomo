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
  #  #Permutations
  #  npermutations=1000,
  #  permutationrun = permutation_test(npermutations, abadult),
  #  permutationstats = stats_permutations(permutationrun, akey, abadult, "akey"),
  #  permutationstats_pey = stats_permutations(permutationrun, pey_coords, abadult, "peycoords")
  #  
   sestan_a = mRNA_sestan(akey),
   sestan_apey = mRNA_sestan(pey_coords),
   #sestan_raw = mRNA_sestan() #Done in cluster
   combinedplot_sestan = plot_sestan_compared(sestan_a, sestan_apey)
)