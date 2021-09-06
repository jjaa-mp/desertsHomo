
# Outlier test
#diffdf = difference_perm_sestan returned object
detect_outlier_tb <- function(stats_dataframe, diffdf, which) {
  test <- stats_dataframe %>% 
    group_by(struct_id, datasource) %>% 
    dplyr::summarize(mean_struct = mean(mean_struct))
  
  test <- dcast(test,  struct_id ~ datasource)
  pdf(paste0(which, "_perm.pdf"))
  plot(test[[2]], test$permutations,
       pch = 16,
       xlab=which,
       ylab="permutations",
       textxy(test[[2]], test$permutations, test$struct_id, offset = -0.4))
  dev.off() #Saving the file
  
  
  
  a <- grubbs.test(diffdf$rest, opposite = FALSE)
  b <- grubbs.test(diffdf$rest, opposite = TRUE)
  
  Matrix::print(diffdf[which.max(diffdf$rest),])
  Matrix::print(a)
  
  Matrix::print(diffdf[which.min(diffdf$rest),])
  Matrix::print(b)
  
} 
