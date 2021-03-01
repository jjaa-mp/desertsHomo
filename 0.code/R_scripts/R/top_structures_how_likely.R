top_structures_how_likely <- function(permutationrun, input, where){
  
  input <- input %>% 
    filter(datasource != "permutations")
  input <- input[order(input$mean_struct,decreasing = TRUE),] # reorder
  #get lower value
  headinput <- input[1,] #highest value stucture
  
  input <- input[order(input$mean_struct,decreasing = FALSE),] # reorder
  tailinput <- input[1,] #lowest value structure
  
  
  #gets bottom structure in permutations:
  test <- permutationrun #permutationrun should come from allperm_order.R
  test <- permrun_order(test)
  expected <- test %>% 
    filter(structure == tailinput$struct_id) 
  expected$probability <- expected$probability/100 #due to ggplot funniness
  expected$position <- stringr::str_remove_all(expected$position, "V")
  expected$position <-  fct_reorder(expected$position,as.integer(expected$position))
  
  bottom <- ggplot(expected, aes(position, probability)) +
    theme_minimal() +
    geom_bar(stat = "identity", fill = "steelblue") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste0("Likelihood of ",expected$structure, " having the lowest expression"),
         subtitle = paste0(where, " relative to control"),
         x="Relative rank per expression mean", 
         y="Probability in control")
  
  #gets bottom structure in permutations:
  expected <- test %>% 
    filter(structure == headinput$struct_id) 
  expected$probability <- expected$probability/100
  expected$position <- stringr::str_remove_all(expected$position, "V")
  expected$position <-  fct_reorder(expected$position,as.integer(expected$position))
  
  
  top <- ggplot(expected, aes(position, probability)) +
    theme_minimal() +
    geom_bar(stat = "identity", fill = "steelblue") +
    scale_y_continuous(labels = scales::percent) +
    labs(title = paste0("Likelihood of ",expected$structure, " having the highest expression"),
         subtitle = paste0(where, " relative to control (random permutations, n = 1000)"),
         x="Relative rank (out of 414, descending)", 
         y="Probability in control")
  
  
  p <- grid.arrange(top,bottom)
  ggsave("top_bottom_str_likelihood.pdf", p, width = 8, height = 8)
  
}


#----
#PLAYGROUND
#----


