test_kruskalwallis <- function(input){
  
  #print("Is the cerebellum specifically different to the rest of tissues in prenatal and child stages?")
  #Extracting relevant data from input
  selected_stages_df <- NULL
  selected_stages_df$prenatal <- input[[1]]
  selected_stages_df$child <- input[[3]]
  
  #shape
  selected_stages_df <- lapply(selected_stages_df, melt)
  selected_stages_df <- ldply(selected_stages_df, data.frame)
  selected_stages_df <-  selected_stages_df %>% 
    dplyr::select(".id", "variable", "value") %>% 
    group_by(variable, .id) %>% 
    dplyr::summarise(mean = mean(value)) 

  #stats
  a <- selected_stages_df %>% 
    filter(.id=="prenatal")
  
  b <- selected_stages_df %>% 
    filter(.id=="child")
  
  #With brain regions as variables
  print("Comparison between child and prenatal stages for all tissues:")
  test <- wilcox.test(a$mean,
                      b$mean,
                      paired=TRUE)
  print(test)
  
  #With cerebellum as variable
  print("Comparison between child and prenatal stages for cerebellum:")
  selected_stages_df <- NULL
  selected_stages_df$prenatal <- input[[1]]
  selected_stages_df$child <- input[[3]]
  
  #shape
  selected_stages_df <- lapply(selected_stages_df, melt)
  selected_stages_df <- ldply(selected_stages_df, data.frame)
  selected_stages_df <-  selected_stages_df %>% 
    dplyr::select(".id", "variable", "value") 
  
  a <- selected_stages_df %>% 
    filter(.id=="prenatal", variable == "CBC_cerebellar cortex")
  
  b <- selected_stages_df %>% 
    filter(.id=="child", variable == "CBC_cerebellar cortex")
  
  #With brain regions as variables
  test <- wilcox.test(a$value,
                      b$value,
                      paired=TRUE)
  print(test)

}