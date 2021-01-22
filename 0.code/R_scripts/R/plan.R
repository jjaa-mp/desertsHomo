# list of functions to call, as stored in .R files

the_plan <-
  drake_plan(
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
  	
  	abakey_data = ab_data_plots(ab, akey, TRUE), # prepares akey data for ab plots
  	abakey_plots = ab_plots(abakey_data, "Akey", " Deserts", "Striatum"),  #generates plots
  	
  	abakeypey_data = ab_data_plots(ab, pey_coords, TRUE), # prepares pey+akey data for ab plots
  	abakepey_plots = ab_plots(abakeypey_data, "AkeyPey", " Deserts and Pey", "Somato - Motor - Parietal - Aud Ctxm"),
  	#generates plots
  	
  	abakeyrac_data = ab_data_plots(ab, rac_coords, TRUE), # prepares rac+akey data for ab plots
  	abakeyrac_plots = ab_plots(abakeyrac_data, "AkeyRac", " Rac and Akey", "Striatum"),  #generates plots
  	
  	rac_coords_noakey = biominput(rac, akey, FALSE), # prepares rac data for ab plots
  	abrac_data = ab_data_plots(ab, rac_coords_noakey, TRUE),
  	abrac_plots = ab_plots(abrac_data, "Rac", " Rac", "Striatum"), #generates plots
  	
  	abakeypey_stats = ab_data_plots(ab, pey_coords, FALSE), 
   abakey_stats = ab_data_plots(ab, akey, FALSE), 
   # prepares akey data for ab plots
  	#Same as abakeypey_data, but without means to calculate stats
  	
   anova_akey = testing_nested_anova(abakey_stats), 
   anova_akeypey = testing_nested_anova(abakeypey_stats), 
   
  	pey_coords_noakey = biominput(pey, akey, FALSE), # prepares pey data for ab plots
  	abpey_data = ab_data_plots(ab, pey_coords_noakey, TRUE),
  	abpey_plots = ab_plots(abpey_data, "Pey", " Pey", "Striatum"),
  	
  	GSEakey = GSE_data(akey, pey_coords, "akey"), # prepares GSE data
  	ttest_result = ttest(GSEakey), 
  	akey_fried = friedman(GSEakey, "akey"), # friedman test for akey
  	#anova_GSEakey = alt_anova(GSEakey), # checking with anova + Tukey if Friedman is alright
   
   
  	GSEpey = GSE_data(akey, pey_coords, "akeypey"),
  	pey_fried = friedman(GSEpey, "pey") # friedman test for pey
   #anova_GSEpey = alt_anova(GSEpey) # checking with anova + Tukey if Friedman is alright
)
