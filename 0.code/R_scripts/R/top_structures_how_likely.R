top_structures_how_likely <- function(permutationrun, input, where){
  
  input <- input %>% 
    filter(datasource != "permutations")
  input <- input[order(input$mean_struct,decreasing = TRUE),] # reorder
  #get lower value
  headinput <- input[1,] #highest value stucture
  
  input <- input[order(input$mean_struct,decreasing = FALSE),] # reorder
  tailinput <- input[1,] #lowest value structure
  
  
  #gets bottom structure in permutations:

  expected <- permutationrun %>% 
    dplyr::filter(structure == tailinput$struct_id) 
    #  %>%  dplyr::mutate(highlight = ifelse(position == "V26", "red", "steelblue")) 
  #currently not used
  expected$probability <- expected$probability/100 #this one due to ggplot funniness
  expected$position <- stringr::str_remove_all(expected$position, "V")
  expected$position <-  fct_reorder(expected$position,as.integer(expected$position))
  
  bottom <- ggplot(expected, aes(position, probability)) +
    theme_minimal() +
    geom_bar(stat = "identity", fill = "steelblue") + 
    #change fill here to 'highlight' column to change
    scale_y_continuous(labels = scales::percent)+
    scale_x_discrete(limits =as.factor(1:26)) +
    labs(title = paste0("Likelihood of ",expected$structure, " having the lowest expression"),
         subtitle = paste0(where, " relative to control"),
         x="Relative rank per expression mean", 
         y="Probability in control")
  
  #gets bottom structure in permutations:
  expected <- permutationrun %>% 
    dplyr::filter(structure == headinput$struct_id) 
    #  %>% dplyr::mutate(highlight = ifelse(position == "V1", "red", "steelblue")) 
  #currently not used
  expected$probability <- expected$probability/100
  expected$position <- stringr::str_remove_all(expected$position, "V")
  expected$position <-  fct_reorder(expected$position,as.integer(expected$position))
  
  
  top <- ggplot(expected, aes(position, probability)) +
    theme_minimal() +
    geom_bar(stat = "identity", fill = "steelblue") + #change fill here to 'highlight' column to change
    scale_y_continuous(labels = scales::percent) +
    scale_x_discrete(limits =as.factor(1:26)) +
    labs(title = paste0("Likelihood of ",expected$structure, " having the highest expression"),
         subtitle = paste0(where, " relative to control"),
         x="Relative rank per expression mean", 
         y="Probability in control") +
    theme(legend.position = "none")
  
  p <- grid.arrange(top,bottom)
  if (where == "Deserts of introgression") {
    name <- c("chen")
    } else {
      name <- c("chenpey")
    }
  
  ggsave(paste0("top_bottom_likelihood", name, ".pdf"), p, width = 8, height = 8)
}



#Notes:
#AA: I tried various t-tests, and everything seems to have a non-random distribution