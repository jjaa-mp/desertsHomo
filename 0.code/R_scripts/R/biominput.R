biominput <- function(bedfile, results) {
  bedfile$V1 <- gsub("chr", "\\1", bedfile$V1) 
  bedfile[2] <- bedfile[2]+1 #Moving to 1-based for Ensembl
  df <- paste(bedfile$V1, bedfile$V2, bedfile$V3, sep = ":") #join
  
  # Using hg19 genome
ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = 'hsapiens_gene_ensembl',
                        host = 'https://grch37.ensembl.org')

  results_bed <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=df,biotype="protein_coding"), mart = ensembl)
  
  both <- results_bed[results_bed$hgnc_symbol %in% results$hgnc_symbol,] #returns overlap
  both <-  both[!(is.na(both$hgnc_symbol) | both$hgnc_symbol==""), ] #Cleaning
  
  return(both)
} 
