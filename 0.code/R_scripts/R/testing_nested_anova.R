testing_nested_anova <- function(input) {

  cerebellum <- lapply(input,'[[',16)
  striatum <- lapply(input,'[[',13)
  
  values <- NULL
  cerebellum <- lapply(cerebellum,'[[',1)
  cerebellum <- reshape2::melt(cerebellum)
  values$cerebellum <- cerebellum$value
  #cerebellum$value because L1 will be duplicated later
  
  striatum <- lapply(striatum,'[[',1)
  striatum <- reshape2::melt(striatum)
  values$striatum <- striatum
  
  values <- as.data.frame(values)
  colnames(values) <- c("cerebellum", "striatum", "stage")
  test <- reshape2::melt(values[1:2])
  test$stage <- values$stage
  
  #stats
  test$stage <- as.factor(test$stage)
  model <- lme(value ~ variable, random=~1|stage,
               data=test,
               method="REML")
  
  anova_results <- anova.lme(model, type="sequential", adjustSigma = FALSE)
  print(anova_results) #Variable p-vale = 0.035
} 
