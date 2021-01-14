abad_5_stages <- function() {
  data("dataset_5_stages")##Loading dataset
  #unique(dataset_5_stages$structure) #Checking structures

  #(Perform enrichment on Genes from Deserts - Skip)
  #input_hyper = data.frame(results$hgnc_symbol, is_candidate=1)
  #res_devel = aba_enrich(input_hyper, dataset='5_stages')

  #Selecting genes ID and structures present in dataset
  id <- unique(dataset_5_stages$ensembl_gene_id)
  st <- unique(dataset_5_stages$structure)
  st_allen <- paste("Allen",st, sep=":") 
  #Expression data for all structures and genes
  ab <- get_expression(structure_ids=st_allen, gene_ids = id, dataset='5_stages')

#Converting data to a dataframe with useful format for later
  list1 = vector(mode="list")
  for (r in 1:length(ab)){
   ab[[r]] <- t(ab[[r]]) #transpose
    list1 <- get_name(colnames(ab[[r]])) #change Allen:XXXX to e.g. M1C_primary motor cortex, etc
    colnames(ab[[r]]) <- list1
    ab[[r]] <- as.data.frame(ab[[r]])
  } 
  return(ab)
}
