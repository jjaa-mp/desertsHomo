allperm_order <- function(permutationrun1, permutationrun2){
  a <- permrun_order(permutationrun1)
  b <- permrun_order(permutationrun2)



  test <- full_join(a, b, c, d, e, f, g, h, i, j)
  test2 <- test %>% 
    group_by(structure, position) %>% 
    dplyr::summarise(probability = mean(probability))
  return(test2)
}