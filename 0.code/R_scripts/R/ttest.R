ttest <- function(dfakey) {
  #dfakey <- dfakey[order(dfakey$max.cluster),]
  #dfakey %>% group_by(organ) %>% mutate(mean1=mean(max.expr), mean2=mean(second.expr))
  meanakey <- dfakey %>% group_by(organ) %>% summarize(Mean = mean(max.expr))
  meanakey[order(meanakey$Mean, decreasing = TRUE),]
  ##Gene expression of genes in Akey - Preliminar test
  result <- pairwise.t.test(dfakey$max.expr, dfakey$organ, p.adjust.method = "BH") 
  print(result)
  return(result)
}