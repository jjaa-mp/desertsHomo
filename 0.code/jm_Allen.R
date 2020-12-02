#Genes within Desert region coordinates
library("biomaRt")

#Using hg19 genome
ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = 'hsapiens_gene_ensembl',
                        host = 'https://grch37.ensembl.org')
#a <- listAttributes(ensembl)
#b <- listFilters(ensembl)


#Akey Deserts coordinates (Ensembl 1-based):
filterlist <- c("1:105400000:120600000", "3:74100000:89300000", "7:106200000:123200000","8:49400000:66500000")

#Select only protein-coding
results=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=filterlist,biotype="protein_coding"), mart = ensembl)

#Extracting gene expression data from Allen Brain Atlas
library(ABAData)
library(ABAEnrichment)
library(dplyr)

#Loading dataset
data("dataset_5_stages")
unique(dataset_5_stages$structure)
unique(dataset_5_stages$age_category)

#Perform enrichment on Genes from Deserts
input_hyper = data.frame(results$hgnc_symbol, is_candidate=1)
res_devel = aba_enrich(input_hyper, dataset='5_stages')

#Selecting genes ID and structures
id <- unique(dataset_5_stages$ensembl_gene_id)
st <- unique(dataset_5_stages$structure)
st_allen <- paste("Allen",st, sep=":") 

#Expression data
ab <- get_expression(structure_ids=st_allen, gene_ids = id, dataset='5_stages')

#Example: Subsetting for Adult gene expression data
ab5 <- t(ab[["age_category_5"]])
  
#Renaming columns
list <- vector(mode = "list")
for (val in 1:length(colnames(ab5))) {
  list[val] <- get_name(colnames(ab5)[val])
}
#print(list)
colnames(ab5) <- list

df2 <- as.data.frame(ab5)

df3 <- df2 %>% select(`DFC_dorsolateral prefrontal cortex`)

f <- df3 %>% filter(df3$`DFC_dorsolateral prefrontal cortex` > quantile(df3$`DFC_dorsolateral prefrontal cortex`, 0.75))
head(rownames(f))

gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = rownames(f), mart= ensembl)

intersect(gene_IDs$hgnc_symbol,results$hgnc_symbol)
