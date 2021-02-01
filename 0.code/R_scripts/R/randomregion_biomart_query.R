randomregion_biomart_query <- function(random) {
  
  ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'https://grch37.ensembl.org')
  
  resultsTemp=getBM(attributes = c("ensembl_gene_id"),
                    filters = c("chromosomal_region","biotype"),
                    values = list(random,biotype="protein_coding"), mart = ensembl)
}