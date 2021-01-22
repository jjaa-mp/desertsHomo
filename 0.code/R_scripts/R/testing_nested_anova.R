testing_nested_anova <- function(input, fulldf) {
  
  if (fulldf == FALSE) {
  print("Using only cerebellum and striatum!")
  #Extracting relevant data from input
  cerebellum <- lapply(input,'[[',16)
  striatum <- lapply(input,'[[',13)
  
  values <- NULL
  cerebellum <- lapply(cerebellum,'[[',1)
  cerebellum <- reshape2::melt(cerebellum)
  values$cerebellum <- cerebellum$value
  #cerebellum$value only because L1 will be duplicated later
  
  striatum <- lapply(striatum,'[[',1)
  striatum <- reshape2::melt(striatum)
  values$striatum <- striatum
  
  values <- as.data.frame(values)
  colnames(values) <- c("cerebellum", "striatum", "stage")
  test <- reshape2::melt(values[1:2])
  test$stage <- values$stage
  #test, the input for the anova, should be tidy now
  }  else if (fulldf == TRUE ) {
    print("Using the whole dataset!")
    df <- input
    df <- lapply(df, melt)
    df <- ldply(df, data.frame)
    stges <- 1:5
    ln <-  length(df$L1)/5
    df$L1 <- rep(stges, each=ln)
    colnames(df)[colnames(df) == "L1"] <- c("stage")
    test <- df
    #test, the input for the anova, should be tidy now
    
  }
  
  #stats
  #With tissues as variables
  test$stage <- as.factor(test$stage)
  model <- lmerTest::lmer(value ~ 1+ variable + (1|stage),
                          data=test)
  anova(model)
  qqnorm(resid(model), main = "Residuals") 
  # Residuals look terrible - model not valid
  
  #Same, but with stages as variables
  # Residuals still look terrible  - model not valid
  test$stage <- as.factor(test$stage)
  model <- lmerTest::lmer(value ~ 1+ stage + (1|variable),
                          data=test)
  anova(model)
  summary(model)
  plot(model)
  qqnorm(resid(model), main = "Residuals") 
}

