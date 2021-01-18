q75 <- function(ab, results) {
	#>q75 to select genes with high expression
	q75 <- vector(mode="list", length = length(ab))
	for (i in 1:length(ab)){
	  for (h in 1:length(names(ab[[i]]))){
	    q75[[i]][[h]] <- ab[[i]][h] %>% 
	      dplyr::filter(ab[[i]][h] > quantile(ab[[i]][[h]], 0.75))
	  }
	}
	
	#Changing row names ENSG to hgnc_symbol (via bioMart); ordering values based on expression
	##ab is a list of 5 elements (i); each element or 'sublist' has 16 dataframes (h)
	##q75[[1]][[3]][[1]] would denote  first sublist (1), third dataframe (3), column 1

	ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = 'hsapiens_gene_ensembl',
                        host = 'https://grch37.ensembl.org')

	for (i in 1:length(ab)){
	  for (h in 1:length(names(ab[[i]]))){
	    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(q75[[i]][[h]]),mart=ensembl)
	    q75[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(q75[[i]][[h]]), G_list$ensembl_gene_id)]
	    q75[[i]][[h]] <-q75[[i]][[h]][order(q75[[i]][[h]][[1]], decreasing = TRUE), ]
	  }
	}
	
	#Intersecting q75 (high expression) with Akey
	akey075 = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
	for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
	  for (h in 1:length(names(ab[[i]]))){
	    akey075[[i]][[h]] <- q75[[i]][[h]][q75[[i]][[h]][[2]] %in% results$hgnc_symbol,] #in Akey
	    akey075[[i]][[h]] <-  akey075[[i]][[h]][!(is.na(akey075[[i]][[h]][[2]]) | akey075[[i]][[h]][[2]]==""), ] #Cleaning
	  }
	}
	
	#Save results - Genes in deserts and pos. selection with high expression
	for (i in 1:length(ab)){
	  for (h in 1:length(akey075[[i]])){
	    if (dim(akey075[[i]][[h]])[1] == 0) next #To skip empty dataframes
	    write.csv(akey075[[i]][[h]], file="output/ABA_Akey_highExprq75.xlsx") #Add age number at the beginning of the sheet name
	  }
	} 
}