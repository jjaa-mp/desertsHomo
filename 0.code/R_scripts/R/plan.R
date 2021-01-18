# list of functions to call, as stored in .R files

the_plan <-
  drake_plan(
  	akey = retrieve_GenesDeserts(),
  	abadult = dt_adults(),
  	
  	pey = load_data("file_dependencies/2020_pey_coords.bed"),
  	pey_coords = biominput(pey, akey, TRUE), # Genes within Akey and Peyregne

  	rac = load_data("file_dependencies/2020_rac_coords.bed"),
  	rac_coords = biominput(rac, akey, TRUE), # Genes within Akey and Racimo

  	abadult_pl = abadult_plots(abadult, 
  	              akey, 
  	              pey_coords),
  	
  	ab = abad_5_stages(),
  	q75(ab, akey),
  	abakey_data = ab_data_plots(ab, akey),
  	abakey_plots = ab_plots(abakey_data, "Akey", " Deserts", "Striatum"),
  	
  	abakeypey_data = ab_data_plots(ab, pey_coords),
  	abakepey_plots = ab_plots(abakeypey_data, "AkeyPey", " Deserts and Pey", "Somato - Motor - Parietal - Aud Ctxm"),
  	
  	abakeyrac_data = ab_data_plots(ab, rac_coords),
  	abakeyrac_plots = ab_plots(abakeyrac_data, "AkeyRac", " Rac and Akey", "Striatum"),
  	
  	rac_coords_noakey = biominput(rac, akey, FALSE),
  	abrac_data = ab_data_plots(ab, rac_coords_noakey),
  	abrac_plots = ab_plots(abrac_data, "Rac", " Rac", "Striatum"),
  	
  	pey_coords_noakey = biominput(pey, akey, FALSE),
  	abpey_data = ab_data_plots(ab, pey_coords_noakey),
  	abpey_plots = ab_plots(abpey_data, "Pey", " Pey", "Striatum"),
  	
  	GSEakey = GSE_data(akey, pey_coords, "akey"),
  	ttest_result = ttest(GSEakey),
  	akey_fried = friedman(GSEakey, "akey"),
  	
  	GSEpey = GSE_data(akey, pey_coords, "akeypey"),
  	pey_fried = friedman(GSEpey, "pey")
  	
)
