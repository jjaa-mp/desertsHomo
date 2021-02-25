top_structures_how_likely <- function(x){
  
  
}

#cgb_cingulum bundle, Right

#----
#PLAYGROUND
#----

loadd(permutationrunABA1)
permrun <- permutationrunABA1

permrun <- melt(permrun)
permrun <- permrun %>% 
  group_by(variable, L1) %>%
  filter(value > 2)
permrun$value <- log2(permrun$value) 

permrun <- permrun %>% 
  group_by(variable, L1) %>% 
  dplyr::summarize(mean_struct = mean(value))

permrun <- permrun %>% 
  group_by(L1) %>% 
  arrange(L1, dplyr::desc(mean_struct), by_group = TRUE)

permrun %>% 
  group_by(L1)