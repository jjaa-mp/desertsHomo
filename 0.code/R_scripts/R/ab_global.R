ab_global <- function(ab, results) {
  ab1 <- ab
  
  ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'https://grch37.ensembl.org')
  
  for (i in 1:length(ab)){
    for (h in 1:length(names(ab[[i]]))){
      ab1[[i]][[h]] <- as.data.frame(ab[[i]][h])
      G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(ab1[[i]][[h]]),mart=ensembl)
      ab1[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(ab1[[i]][[h]]), G_list$ensembl_gene_id)]
      ab1[[i]][[h]] <-ab1[[i]][[h]][order(ab1[[i]][[h]][[1]], decreasing = TRUE), ]
    }
  }
  
  
  aba_akey = vector(mode="list", length = length(ab1))
  for (i in 1:length(ab1)){
    for (h in 1:length(names(ab1[[i]]))){
      aba_akey[[i]][[h]] <- ab1[[i]][[h]][ab1[[i]][[h]][[2]] %in% results$hgnc_symbol,] #in Akey
      aba_akey[[i]][[h]] <-  aba_akey[[i]][[h]][!(is.na(aba_akey[[i]][[h]][[2]]) | aba_akey[[i]][[h]][[2]]==""), ] #Cleaning
    }
  }
  
  #Reporting mean
  mean1 = vector(mode="list", length = length(ab1))
  for (i in 1:length(ab1)){
    for (h in 1:length(names(ab1[[i]]))){
      mean1[[i]][[h]] <- mean(aba_akey[[i]][[h]][[1]])
      mean1[[i]][[h]] <- as.data.frame(mean1[[i]][[h]])
      names(mean1[[i]][[h]]) <- paste(names(aba_akey[[i]][[h]][1]), sep='_')
    }
  }
  
  prenatal <- do.call("rbind", as.data.frame(mean1[[1]]))
  colnames(prenatal) <- "prenatal"
  infant <- do.call("rbind", as.data.frame(mean1[[2]]))
  colnames(infant) <- "infant"
  child <- do.call("rbind", as.data.frame(mean1[[3]]))
  colnames(child) <- "child"
  adolescent <- do.call("rbind", as.data.frame(mean1[[4]]))
  colnames(adolescent) <- "adolescent"
  adult <- do.call("rbind", as.data.frame(mean1[[5]]))
  colnames(adult) <- "adult"
  
  final_merge <- as.data.frame(cbind(prenatal, infant, child, adolescent, adult))
  #Save results - need to be updated
  library(xlsx)
  for (i in 1:length(colnames(final_merge))){
    write.csv(arrange(final_merge[i], desc(final_merge[i])), file="output/csv")
  }
  
  #Plot
  
  #For plot
  final_merge <- tibble::rownames_to_column(final_merge, "Structure")
  final_merge[[1]] <- sapply(strsplit(final_merge[[1]], "_"), "[", 1)
  
  a<-ggparcoord(final_merge,
                columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts")+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("mean  expression")
  b <- ggparcoord(final_merge,
                  columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("mean expression")
  c<-ggparcoord(final_merge,
                columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("mean expression")
  a1<-arrangeGrob(a, left=textGrob("A"))
  b1<-arrangeGrob(b, left =textGrob("B"))
  c1<-arrangeGrob(c, left=textGrob("C"))
  grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)
  
  #pdf("ABA_GenesAkey.pdf", paper="a4")
  #grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
  #dev.off()
  
  pl1 <-grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
  ggsave(file="output/ABA_GenesAkey.pdf", pl1, width = 11.69, height = 8.27, units = "in") 
}