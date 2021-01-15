ab_plots <- function(final_merge, out, title, struct) {

	 #Plot
  final_merge <- tibble::rownames_to_column(final_merge, "Structure")
  final_merge[[1]] <- sapply(strsplit(final_merge[[1]], "_"), "[", 1)
  
  a<-ggparcoord(final_merge,
                columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title=paste0("Genes in", title))+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("mean  expression")
  b <- ggparcoord(final_merge,
                  columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title=struct)+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("mean expression")
  c<-ggparcoord(final_merge,
                columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("mean expression")
  a1<-arrangeGrob(a, left=textGrob("A"))
  b1<-arrangeGrob(b, left =textGrob("B"))
  c1<-arrangeGrob(c, left=textGrob("C"))
  grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)
  
  #pdf("ABA_GenesAkey.pdf", paper="a4")
  #grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
  #dev.off()
  
  if (struct == "Striatum"){
  pl1 <-grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
  ggsave(file=paste0("output/ABA_Genes",out,".pdf"), pl1, width = 11.69, height = 8.27, units = "in") 
  #ggsave(file=paste0("../../2.plots/ABAData_AkeyPeyRac/ABA_Genes",out,".pdf"), pl1, width = 11.69, height = 8.27, units = "in") 
  } else {
    d<-ggparcoord(final_merge,
                  columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Somato - Motor - Parietal - Aud Ctx")+scale_color_manual(values = c( "#00FF00", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#00FF00", "#ABABAB",  "#00FF00", "#ABABAB", "#ABABAB","#ABABAB", "#00FF00", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("mean expression")
    d1<-arrangeGrob(d, left=textGrob("D"))
    
    pl2 <- grid.arrange(a1, c1, d1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
    ggsave(file="output/ABA_GenesAkeyPey.pdf", pl2, width = 11.69, height = 8.27, units = "in")
    
  }

}
