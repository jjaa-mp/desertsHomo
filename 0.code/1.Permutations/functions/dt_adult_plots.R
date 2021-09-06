dt_adults <- function() {
  data("dataset_adult") #from ABAD package
  # Getting data
  id.adult <- unique(dataset_adult$ensembl_gene_id)
  st.adult <- unique(dataset_adult$structure)
  st.adult <- paste("Allen",st.adult, sep=":") 
  abadult <- get_expression(structure_ids=st.adult, gene_ids = id.adult, dataset='adult') 
  abadult <- t(abadult)
  listadult <- vector(mode="list")
  listadult <- get_name(colnames(abadult))
  colnames(abadult) <- listadult
  abadult <- as.data.frame(abadult)

  # Using hg19 genome
  ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'https://grch37.ensembl.org')
  
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",  "hgnc_symbol"),values=rownames(abadult),mart=ensembl)

  abadult['gene_name'] <-  G_list$hgnc_symbol[match(rownames(abadult), G_list$ensembl_gene_id)]
  nrow(abadult) # 15698

  abadult <- abadult[!(is.na(abadult$gene_name) | abadult$gene_name==""), ]
  nrow(abadult) # 15688

  row.names(abadult) <- G_list$hgnc_symbol[match(rownames(abadult), G_list$ensembl_gene_id)]
  return(abadult)
}
