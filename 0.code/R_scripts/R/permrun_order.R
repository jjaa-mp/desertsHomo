permrun_order <- function(orig_permrun){
  
  permrun <- orig_permrun # to keep track of n.er of permutations
  permrun <- melt(permrun)
  permrun <- permrun %>% 
    group_by(variable, L1) %>%
    filter(value > 2)
  permrun$value <- log2(permrun$value) 
  
  
  permrun$variable <- stringr::str_remove_all(permrun$variable, ".*[.]")
  permrun <- permrun %>% 
    group_by(variable, L1) %>% 
    dplyr::summarize(mean_struct = mean(value))
  
  
  
  permrun <- permrun %>% 
    group_by(L1) %>% 
    arrange(L1, dplyr::desc(mean_struct), by_group = TRUE)
  
  dforder <- NULL
  s <- 1 #start row
  e <-26 #end row
  colnum <- 1
  while (e < nrow(permrun)+1) {
    dforder[[colnum]] <- permrun$variable[s:e]
    s <- s+26
    e <- e+26
    colnum <- colnum+1
  }
  
  dforder <- as.data.frame(dforder)
  colnames(dforder) <- c(1:ncol(dforder))
  
  dforder <- t(dforder) 
  # dforder = each row is a permutation run (there are 100, minus empty ones, usually ~92)
  # each column, the rank in expression of that structure (there are 414)
  
  test <- as.data.frame(dforder)
  test <- lapply(test, table) 
  test <- melt(test)
  test$value <- test$value/length(orig_permrun) #percentage of times in a given position
  #test now sets the controlled, expected probabilities 
  
  colnames(test) <- c("structure", "probability", "position")
  return(test)
}