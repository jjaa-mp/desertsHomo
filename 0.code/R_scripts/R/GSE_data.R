GSE_data <- function (results, relevantregion, out) {
  df <- read.csv(file = "file_dependencies/GSE156793_S8_DE_gene_cells.csv.gz")
  df1 <- df[which(df$organ=='Cerebellum'), ]
  brain <- df[which(df$organ=='Cerebrum'), ]
  write.csv(df1, file="output/CellAtlas_GSE156793_cerebellum.csv", row.names = FALSE)
  
  df1$gene_short_name <- gsub("\\'", "", df1$gene_short_name)
  df1subset <- df1[df1$gene_short_name %in% results$hgnc_symbol,]
  df1subset <- df1subset[order(df1subset$max.cluster),]
  
  write.csv(df1subset, file="output/CellAtlas_GSE156793_inAkey.csv", row.names = FALSE)
  df$gene_short_name <- gsub("\\'", "", df$gene_short_name)
  
  if (out == "akey"){
  #Whole dataset - AKEY
  dfakey <- df[df$gene_short_name %in% results$hgnc_symbol,]
  return(dfakey)
} else if (out == "akeypey"){
    
    df1subsetboth <- df1[df1$gene_short_name %in% relevantregion$hgnc_symbol,]
    df1subsetboth <- df1subsetboth[order(df1subsetboth$max.cluster),]
    df1subset <- df1subset[order(df1subset$max.cluster),]
    write.csv(df1subsetboth, file="output/CellAtlas_GSE156793_inAkeyPey.csv", row.names = FALSE)
    #Whole dataset - AKEY PEY
    df1subsetboth <- df[df$gene_short_name %in% relevantregion$hgnc_symbol,]
    #meanakeypey <- df1subsetboth %>% group_by(organ) %>% summarize(Mean = mean(max.expr), .groups = 'drop')
    #meanakeypey[order(meanakeypey$Mean, decreasing = TRUE),]
    return(df1subsetboth)
    
  } else{
    print("Something went wrong")
  }


  #Raw:
  #df  %>% group_by(organ) %>% summarize(Mean = mean(max.expr), .groups = 'drop')
  
} 
