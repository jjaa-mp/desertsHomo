nemenyi_posthoc <- function(df, outputname) {
  
  if (outputname == "chen") {
    filename <- "chen"
  } else if (outputname == "chenpey") {
    filename <- "chen_pey"
  }
  
  # Currently drake complains about write.csv - not sure why still. 
  # Recommended to run manually with:
  # loadd(abakey_stats) OR loadd(abakeypey_stats)
  # df <- abakey_stats / abakeypey_stats
  # Then the following for loop:
    for (stage in 1:5){
      test <- df[[stage]]
      test <- melt(test)
      ph <- frdAllPairsNemenyiTest(value ~ variable | gene_name, data = test)
      #write.csv(ph$p.value, file = paste0("output/", outputname, stage, ".csv", ))
      
      test2 <- ph$p.value
      colnames(test2) <- stringr::str_remove_all(colnames(test2), regex("_.*"))
      rownames(test2) <- stringr::str_remove_all(rownames(test2), regex("_.*"))
      test2 <- melt(test2)
      
      plot <- ggplot(data = test2, aes(Var2, Var1, fill = value))+
        theme_minimal() +
        theme(panel.grid = element_blank()) +
        geom_raster()+
        scale_fill_gradient2(low = "#F0E442", high = "#0072B2", mid = "#CC79A7", 
                             midpoint = 0.5, limit = c(0,1), space = "Lab", 
                             name="p-value", na.value = 'white') +
        xlab("")+
        ylab("")+
        theme(axis.text.x = element_text(angle = -90,hjust = 0))+
        scale_y_discrete(position = "left") + 
        coord_fixed()
      ggsave(file=paste0("../../2.plots/ABAData_AkeyPeyRac/posthoc_nemenyi/", outputname, "_", stage, ".pdf"),plot,width = 12, height = 12, units = "in")
    }
}
