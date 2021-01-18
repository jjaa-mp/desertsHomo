friedman <- function(df, title) {
  r <- data.frame(df$organ, df$gene_short_name, df$max.expr)
  names(r) <- c("x", "y", "z")
  rr <- pivot_wider(r, names_from = x, values_from = z)
  rr3 <- as.matrix(rr)
  rr3 <-rr3[,-1]
  friedman.test(rr3)
  
  
  list1 <- vector(mode = "list")
  for (i in 1:length(colnames(rr3))){
    list1[[i]] <- as.numeric(c(rr3[,i]))
  }
  res_F <- DunnettTest(list1, control = c(1:15)) #CAN ALSO TRY SNK Test
  organs <- colnames(rr3)
  list_dfs <- vector(mode="list")
  for (i in 1:length(res_F)){
    colnames(res_F[[i]]) <- c("diff", "lwr.ci", "upr.ci", colnames(rr3)[i]) #control group: colnames(rr3)[i]
    rownames(res_F[[i]]) <- organs[-i] #comparison groups in rownames: all except i
    list_dfs[[i]] <- as.data.frame(res_F[[i]])
    list_dfs[[i]] <- list_dfs[[i]][4]
    list_dfs[[i]] <-  cbind(list_dfs[[i]], group=rownames(list_dfs[[i]]))
    list_dfs[[i]] <- list_dfs[[i]][c(2,1)] #ordering columns for merging on index 'group'
  }
  
  pp <-reduce(list_dfs, full_join, by="group") #Generating matrix for plot
  pp[is.na(pp)] <- 2 #replacing null values
  rownames(pp) <- pp[,1]
  pp[,1] <- NULL
  pp <- pp[c(2:15,1)] #Reordering columns for triangular matrix
  
  upper_tri <- get_upper_tri(as.matrix(pp))
  melted <- melt(upper_tri, na.rm = TRUE)
  
  pl_fr <- ggplot(data = melted, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "#F0E442", high = "#0072B2", mid = "#CC79A7", 
                         midpoint = 0.5, limit = c(0,1), space = "Lab", 
                         name="p-value") +
    theme_classic()+xlab("")+ylab("")+ 
    theme(axis.text.x = element_text(angle = 45,hjust = 0))+scale_y_discrete(position = "right")+
    coord_fixed()+ coord_flip()
  pl_fr
  ggsave(file=paste0("output/CellAtlas_MeanExpr", title, ", _heatmap.pdf"), pl_fr,width = 11.69, height = 8.27, units = "in")
  
} 
