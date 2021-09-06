plotting_topbottom <- function (full_data, whichone, out) {
if (whichone == "top") {
  plot <- full_data %>% 
    filter(datasource != "permutations") %>% 
    unique() %>% 
    group_by(struct_id) %>% 
    #derive minimum value per stucture
    slice_max(order_by = mean_struct, n = 1) %>% 
    ungroup(struct_id) %>% 
    slice_max(order_by = mean_struct, n = 20)
} else if (whichone == "bottom") {
  #Some structures are repeated, so this is going to be a bit more elaborated
    plot <- full_data %>% 
      filter(datasource != "permutations") %>% 
      unique() %>%
      group_by(struct_id) %>% 
      #derive minimum value per stucture
      slice_min(order_by = mean_struct, n = 1) %>% 
      ungroup(struct_id) %>% 
      slice_min(order_by = mean_struct, n = 20)
}
  
other <- full_data %>% 
  filter(struct_id %in% plot$struct_id) %>% 
  filter(datasource == "permutations") %>% 
  group_by(struct_id) %>% 
  dplyr::summarize(mean_struct = mean(mean_struct)) %>% 
  dplyr::mutate(datasource = "permutations")

new <- rbind(plot,other)
p <- ggplot(new, aes(x = mean_struct, y = struct_id, colour = datasource)) +
  theme_minimal() +
  geom_point() +
  labs(x = "Log2 of mean expression", y = "Structure", colour = "Source of data" )

if (out == "chen") {
  ggsave(paste0(whichone, "20chen.pdf"), p, width = 8, height = 4)
} else if (out == "chenpey") {
  ggsave(paste0(whichone, "20chenpey.pdf"), p, width = 8, height = 4)
}

} 