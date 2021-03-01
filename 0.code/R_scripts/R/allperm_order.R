allperm_order <- function(orderperm1, orderperm2,
                          orderperm3, orderperm4,
                          orderperm5, orderperm6,
                          orderperm7, orderperm8,
                          orderperm9, orderperm10) {

  test <- rbind(orderperm1, orderperm2,
                        orderperm3, orderperm4,
                        orderperm5, orderperm6,
                        orderperm7, orderperm8,
                        orderperm9, orderperm10)
  test <- test %>% 
    group_by(structure, position) %>% 
    dplyr::summarise(probability = mean(probability))
  return(test)
}