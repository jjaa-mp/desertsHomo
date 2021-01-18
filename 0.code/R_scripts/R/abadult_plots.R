abadult_plots <- function (abadult, results, pey_coords){
  #In Akey:
  abadultAkey  <- abadult[rownames(abadult) %in% results$hgnc_symbol,]
  abadultAkeyPey  <- abadult[rownames(abadult) %in% pey_coords$hgnc_symbol,]
  
  new <- as.data.frame(colMeans(abadultAkey[sapply(abadultAkey, is.numeric)]))
  names(new) <- "mean_expression"
  
  newboth <- as.data.frame(colMeans(abadultAkeyPey[sapply(abadultAkeyPey, is.numeric)]))
  names(newboth) <- "mean_expression"
  newboth <- arrange(newboth, desc(newboth$mean_expression))
  newboth1 <- newboth %>% slice(head(row_number(), 20)) #Top 20 structures
  p<-ggplot(newboth1, aes(x=rownames(newboth1), y=newboth1$mean_expression)) + 
    geom_dotplot(binaxis='y', stackdir='center', fill="#D55E00")+theme(legend.position = "none")+labs(title="",x="", y = "Mean expression (top 20)")+coord_flip()
  
  newboth2 <- newboth %>% slice(tail(row_number(), 20)) #Bottom 20 structures
  p2<-ggplot(newboth2, aes(x=rownames(newboth2), y=newboth2$mean_expression)) + 
    geom_dotplot(binaxis='y', stackdir='center', fill="#0072B2")+theme(legend.position = "none")+labs(title="",x="", y = "Mean expression (bottom 20)") + coord_flip()+scale_x_discrete(position = "top")
  print("plots ready!")
  ggsave(file="output/ABA_414_GenesAkeyPey_top20.pdf", p, width = 11.69, height = 8.27, units = "in")
  ggsave(file="output/output_ABA_414_GenesAkeyPey_bottom20.pdf", p2, width = 11.69, height = 8.27, units = "in") 
} 
