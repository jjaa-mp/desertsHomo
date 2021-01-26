test_nonparam_anova <- function(df, output) {
  df <- lapply(df, melt)
  df <- ldply(df, data.frame)
  
  #stages
  stages <- 1:5
  ln <-  length(df$L1)/5
  df$L1 <- rep(stages, each=ln)
  colnames(df)[colnames(df) == "L1"] <- c("stage")
  df$stage <- str_replace_all(df$stage, "1", "Prenatal")
  df$stage <- str_replace_all(df$stage, "2", "Infant")
  df$stage <- str_replace_all(df$stage, "3", "Child")
  df$stage <- str_replace_all(df$stage, "4", "Adolescent")
  df$stage <- str_replace_all(df$stage, "5", "Adult")

  df$stage <- as.factor(df$stage)
  model <- art(value ~ variable*stage, data=df)
  summary(model)
  print(anova(model))
  
  # Constrasts tests 
  print("Contrast of main effect: Stages only")
  res <- contrast(emmeans(artlm(model, "stage"), ~ stage), method="pairwise")
  print(summary(res))
  write.csv(res, file = paste0("output/anova/", output, "_stages.csv"))
  
  print("Contrast of main effect: organs only")
  res <- contrast(emmeans(artlm(model, "variable"), ~ variable), method="pairwise")
  print(summary(res))
  write.csv(res, file = paste0("output/anova/", output, "_organs.csv"))
  
  print("Contrast of main effects: organs and stages")
  res <- contrast(emmeans(artlm(model, "variable:stage"), ~ variable:stage), method="pairwise", interaction=TRUE)
  write.csv(res, file = paste0("output/anova/", output, "_stages_organs.csv"))
}
