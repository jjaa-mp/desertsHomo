retrieve_GenesDeserts <- function() {
	# Import files (with start 0-based)
	pey_coords <- read.delim("file_dependencies/2020_pey_coords.bed", header=FALSE) 
	rac_coords <- read.delim("file_dependencies/2020_rac_coords.bed", header=FALSE)

	# Using hg19 genome
	ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = 'hsapiens_gene_ensembl',
                        host = 'https://grch37.ensembl.org')

	# Akey Deserts coordinates (Ensembl 1-based):
	filterlist <- c("1:105400000:120600000", "3:74100000:89300000", "7:106200000:123200000","8:49400000:66500000")

	# Select only protein-coding from Akey
	results <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=filterlist,biotype="protein_coding"), mart = ensembl)

	#write.csv(results, file="./output/2020_Genes_in_Deserts.csv") 
	return(results)
}


# Alternative: 800 genes (protein coding genes plus other genes)
#results_alternative <- getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
#              filters = c("chromosomal_region"),
#              values = list(chromosomal_region=filterlist), mart = ensembl)