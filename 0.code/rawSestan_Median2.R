mRNAseqData=read.table("mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
modsb1= mRNAseqData %>% 
  separate(Geneid,c("EnsemblID","Genename"),extra="merge")
modsb1$EnsemblID<-NULL
rawSestan <- modsb1
rawSestan=t(rawSestan)
#As dataframe
rawSestan=as.data.frame(rawSestan)
colnames(rawSestan) <- as.matrix(unlist(rawSestan[1,]))
rawSestan <- rawSestan[-1, ]
rawSestan <- cbind(info = rownames(rawSestan), rawSestan)
rownames(rawSestan) <- 1:nrow(rawSestan)
#duplicated columns - Remove duplicates if needed
#colnames(rawSestan)[duplicated(colnames(rawSestan))] #0
rawSestan <- rawSestan[, !duplicated(colnames(rawSestan))]
rawSestan=rawSestan %>% 
  separate(info, c("Braincode","Regioncode"))
lograwsestan1 <- rawSestan
cols.num <- colnames(lograwsestan1[3:ncol(lograwsestan1)])
lograwsestan1[,cols.num] <- lapply(lograwsestan1[cols.num],as.character)
lograwsestan1[,cols.num] <- lapply(lograwsestan1[cols.num],as.numeric)
lograwsestan1[3:ncol(lograwsestan1)] <- log2(lograwsestan1[3:ncol(lograwsestan1)]+1)
metadatamRNAseq=read_xlsx("mRNA-seq_QC.xlsx",skip = 3)
modMetadatamRNAseq=na.omit(metadatamRNAseq)
modMetadatamRNAseq=modMetadatamRNAseq %>% select(1:3)
modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)
finalrawSestan=merge(modMetadatamRNAseq,lograwsestan1,by=c("Braincode", "Regioncode"))
finalraw2 <- finalrawSestan
rownames(finalraw2)<- do.call(paste,c(finalraw2[c("Braincode","Regioncode", "Window")],sep="_"))
finalraw2[(1:3)] <- NULL
finalraw2<-finalraw2 %>% select(which(colMedians(as.matrix(finalraw2))>2))
finalraw2 <- rownames_to_column(finalraw2) 
finalraw2 <- finalraw2 %>% 
  separate(rowname, c("Braincode","Regioncode", "Window"))
finalraw2$Braincode <- NULL
finalraw2$Label<- paste(finalraw2$Regioncode, finalraw2$Window, sep="_")
finalraw2$Regioncode <- NULL
finalraw2$Window <- NULL
finalraw2 <- as_tibble(finalraw2)
finalraw2_w_mean <- finalraw2 %>%
  group_by(Label) %>%
  dplyr::summarise_all(median, na.rm=TRUE)
dfrawSestan2 <- data.frame(Structure=finalraw2_w_mean$Label, Means=rowMeans(finalraw2_w_mean[,-1]))
dfrawSestan2 <- dfrawSestan2 %>% 
  separate(Structure, c("Structure","Window"))
df_raw2 <- pivot_wider(dfrawSestan2, names_from = Window, values_from = Means)
df_raw2[,10] <- NULL
df_raw2 <- df_raw2[complete.cases(df_raw2), ]
colnames(df_raw2) <- c("Structure", "Fetal1", "Fetal2", "Fetal3", "Birth/Ifan", "Infan/Childh", "Childh", "Adolescence", "Adulth")
#PLOT
levels(colnames(df_raw2)) <- c("Structure", "Fetal1", "Fetal2", "Fetal3", "Birth/Ifancy", "Infancy/Childh", "Childh", "Adolescence", "Adulth")
#write.csv(df_raw2, file="median_filtered_rawSestan.csv", row.names = FALSE)
a<- ggparcoord(df_raw2,
               columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts- VFC (green) & AMY (black)")+scale_color_manual(values = c( "#ABABAB", "#000000", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#238b45"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
b <- ggparcoord(df_raw2,
                columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
c<-ggparcoord(df_raw2,
              columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("expression")
a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
dfraw2_pl <-grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
#ggsave(file="filtered_rawSestan_log2_median.pdf", dfraw2_pl, width = 11.69, height = 8.27, units = "in")

