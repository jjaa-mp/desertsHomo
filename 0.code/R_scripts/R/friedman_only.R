friedman_only <- function(df){
  df <- melt(df)
  
  r <- data.frame(df$variable, df$gene_name, df$value)
  names(r) <- c("x", "y", "z")
  r <- pivot_wider(r, names_from = x, values_from = z)
  r <- as.matrix(r)
  r <- r[,-1]
  friedman.test(r)
  print(paste0("Bonferroni correction: ", (0.05)/5))
}
