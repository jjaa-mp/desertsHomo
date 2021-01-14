# list of functions to call, as stored in .R files

the_plan <-
  drake_plan(
  	results = retrieve_GenesDeserts(),
  	abadult = dt_adults(),
  	
  	pey = load_data("file_dependencies/2020_pey_coords.bed"),
  	pey_coords = biominput(pey, results), # Genes within Akey and Peyregne

  	rac = load_data("file_dependencies/2020_rac_coords.bed"),
  	rac_coords = biominput(rac, results), # Genes within Akey and Racimo

  	abadult_pl = abadult_plots(abadult, 
  	              results, 
  	              pey_coords),
  	
  	ab = abad_5_stages(),
  	q75(ab, results),
  	ab_global(ab, results)
)
