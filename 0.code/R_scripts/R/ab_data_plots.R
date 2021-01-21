ab_data_plots <- function(ab, where, means) {
  ab1 <- ab
  
  ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'https://grch37.ensembl.org')
  
  for (i in 1:length(ab)){
    for (h in 1:length(names(ab[[i]]))){
      ab1[[i]][[h]] <- as.data.frame(ab[[i]][h])
      G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(ab1[[i]][[h]]),mart=ensembl)
      ab1[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(ab1[[i]][[h]]), G_list$ensembl_gene_id)]
      ab1[[i]][[h]] <-ab1[[i]][[h]][order(ab1[[i]][[h]][[1]], decreasing = TRUE), ]
    }
  }
  
  
  aba_akey = vector(mode="list", length = length(ab1))
  for (i in 1:length(ab1)){
    for (h in 1:length(names(ab1[[i]]))){
      aba_akey[[i]][[h]] <- ab1[[i]][[h]][ab1[[i]][[h]][[2]] %in% where$hgnc_symbol,] 
      aba_akey[[i]][[h]] <-  aba_akey[[i]][[h]][!(is.na(aba_akey[[i]][[h]][[2]]) | aba_akey[[i]][[h]][[2]]==""), ] #Cleaning
    }
  }
  
  #Reporting mean
  mean1 = vector(mode="list", length = length(ab1))
  for (i in 1:length(ab1)){
    for (h in 1:length(names(ab1[[i]]))){
      mean1[[i]][[h]] <- mean(aba_akey[[i]][[h]][[1]])
      mean1[[i]][[h]] <- as.data.frame(mean1[[i]][[h]])
      names(mean1[[i]][[h]]) <- paste(names(aba_akey[[i]][[h]][1]), sep='_')
    }
  }
  
  if (means == TRUE){
  prenatal <- do.call("rbind", as.data.frame(mean1[[1]]))
  colnames(prenatal) <- "prenatal"
  infant <- do.call("rbind", as.data.frame(mean1[[2]]))
  colnames(infant) <- "infant"
  child <- do.call("rbind", as.data.frame(mean1[[3]]))
  colnames(child) <- "child"
  adolescent <- do.call("rbind", as.data.frame(mean1[[4]]))
  colnames(adolescent) <- "adolescent"
  adult <- do.call("rbind", as.data.frame(mean1[[5]]))
  colnames(adult) <- "adult"
  final_merge <- as.data.frame(cbind(prenatal, infant, child, adolescent, adult))
  return(final_merge)
  }
  else if (means == FALSE){
    return(aba_akey)
  }
 #Save results - need to be updated
  #for (i in 1:length(colnames(final_merge))){
  #  write.xlsx(arrange(final_merge[i], desc(final_merge[i])), file=paste0("ABA_Genes", where, ".xlsx"), sheetName=names(final_merge[i]), append = TRUE)
  #}
  

}
