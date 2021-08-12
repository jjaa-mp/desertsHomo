#Internal function for permutation_test.R
filter_genexpr_random <- function(region, abadult) {
  subset.data <- abadult %>%
    dplyr::filter(gene_name %in% region$hgnc_symbol) %>% 
    dplyr::select(-c(gene_name))
  return(subset.data) 
}
