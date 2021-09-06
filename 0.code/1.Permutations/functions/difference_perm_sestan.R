
# Plot 1
difference_perm_sestan <- function(stats_perm_df, where, which, height){
  
  stats_perm_df <- stats_perm_df %>% 
    group_by(struct_id, datasource) %>% 
    dplyr::summarize(mean_struct = mean(mean_struct))
  
  stats_perm_df$struct_id <- as.factor(stats_perm_df$struct_id)
  
  #stacked
  ggplot(stats_perm_df, aes (x= mean_struct, y = struct_id, fill = datasource )) +
    theme_minimal() + 
    geom_histogram(stat = "identity")
  
  diff <- stats_perm_df %>% 
    group_by(struct_id) %>% 
    dplyr::mutate(rest = mean_struct[datasource != 'permutations'] - mean_struct[datasource == 'permutations'] ) %>% 
    filter(datasource == 'permutations')
  
  #reorder
  diff$struct_id <- factor(diff$struct_id,                                    # Factor levels in decreasing order
                           levels = diff$struct_id[order(diff$rest, decreasing = FALSE)])
  
  a <-  ggplot(diff, aes (x= rest, y = struct_id )) +
    theme_minimal() + 
    labs(x = c("Log2 mean difference"), y = c("Structure")) +
    geom_histogram(stat = "identity",  color="steelblue", fill="steelblue") 
  
  ggsave(paste0("difference_perm_", where, which, ".pdf"), width = 8, height = height)
  
  return(diff)
} 
