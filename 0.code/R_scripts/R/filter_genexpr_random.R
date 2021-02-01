#Internal function for permutation_test.R
filter_genexpr_random <- function(region) {
  subset.data <- abadult %>%
    dplyr::filter(gene_name %in% region$hgnc_symbol) %>% 
    dplyr::select(-c(gene_name))
  subset.data <- subset.data %>%
    filter(. > quantile(0.10))
  return(subset.data) 
}
