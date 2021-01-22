alt_anova <- function(df) {
  r <- data.frame(df$organ, df$gene_short_name, df$max.expr)
  r$df.organ <- as.factor(r$df.organ)
  model <- lme(df.max.expr ~ 1+ df.organ, random=~1|df.gene_short_name,
               data=r,
               method="REML")
  
  anova_results <- anova.lme(model, type="sequential", adjustSigma = FALSE)
  print(anova_results) #Variable p-vale = 0.035
  
  posthoc <- glht(model,
                 linfct = mcp(df.organ="Tukey"), test = adjusted(type = "bonferroni"))
  print(summary(posthoc))
}


