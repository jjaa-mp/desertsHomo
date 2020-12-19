```{r}
#Packages
library(biomaRt)
library(ABAData)
library(ABAEnrichment)
library(dplyr)
library(GGally)
library(viridis)
library(ggpubr)
library(grid)
library(gridExtra)
library(lattice)
library(xlsx)
library(tidyr)
library(DescTools)
library(reshape2)
```

#Extracting gene names from coordinate regions via bioMart
```{r}
#Genes within Desert region coordinates
##Using hg19 genome
ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = 'hsapiens_gene_ensembl',
                        host = 'https://grch37.ensembl.org')
##Akey Deserts coordinates (Ensembl 1-based):
filterlist <- c("1:105400000:120600000", "3:74100000:89300000", "7:106200000:123200000","8:49400000:66500000")
##Select only protein-coding from Akey
results=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=filterlist,biotype="protein_coding"), mart = ensembl)

#ALTERNATIVE: 800 genes (protein coding genes plus other genes)
##results_alternative =getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
              #filters = c("chromosomal_region"),
              #values = list(chromosomal_region=filterlist), mart = ensembl)
##Save:
#write.csv(results, file="2020_Genes_in_Deserts.csv")

##Pey coordinates (Ensembl 1-based):
pey_coords <- read.delim("~/2020_pey_coords.bed", header=FALSE) #File with start 0-based
###Preparing bed file for input in bioMart
pey_coords$V1 <- gsub("chr", "\\1", pey_coords$V1)
pey_coords[2] <- pey_coords[2]+1 #Moving to 1-based for Ensembl
df <- paste(pey_coords$V1, pey_coords$V2, pey_coords$V3, sep = ":")
results_pey=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=df,biotype="protein_coding"), mart = ensembl)
###Genes within Akey and Pey - RESULT
both <- results_pey[results_pey$hgnc_symbol %in% results$hgnc_symbol,]
both <-  both[!(is.na(both$hgnc_symbol) | both$hgnc_symbol==""), ] #Cleaning

###Rac coordinates #MOVE TO ENSEMBL 1 BASED
rac_coords <- read.delim("~/2020_rac_coords.bed", header=FALSE)
rac_coords$V1 <- gsub("chr", "\\1", rac_coords$V1)
rac_coords[2] <- rac_coords[2]+1

df1 <- paste(rac_coords$V1, rac_coords$V2, rac_coords$V3, sep = ":")
results_rac=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=df1,biotype="protein_coding"), mart = ensembl)

###Genes within Akey and Rac - RESULT
racAkey <- results_rac[results_rac$hgnc_symbol %in% results$hgnc_symbol,]
racAkey <- racAkey[!(is.na(racAkey$hgnc_symbol) | racAkey$hgnc_symbol==""), ]
```

#Extracting gene expression data from Allen Brain Atlas - ADULT
```{R}
data("dataset_adult")
id.adult <- unique(dataset_adult$ensembl_gene_id)
st.adult <- unique(dataset_adult$structure)
st.adult <- paste("Allen",st.adult, sep=":") 
abadult <- get_expression(structure_ids=st.adult, gene_ids = id.adult, dataset='adult') 
abadult <- t(abadult)
listadult = vector(mode="list")
listadult <- get_name(colnames(abadult))
colnames(abadult) <- listadult
abadult <- as.data.frame(abadult)

G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(abadult),mart=ensembl)
abadult['gene_name'] <-  G_list$hgnc_symbol[match(rownames(abadult), G_list$ensembl_gene_id)]
nrow(abadult) # 15698
abadult <- abadult[!(is.na(abadult$gene_name) | abadult$gene_name==""), ]
nrow(abadult) # 15688
row.names(abadult) <- G_list$hgnc_symbol[match(rownames(abadult), G_list$ensembl_gene_id)]
#In Akey:
abadultAkey  <- abadult[rownames(abadult) %in% results$hgnc_symbol,]
abadultAkeyPey  <- abadult[rownames(abadult) %in% both$hgnc_symbol,]

new <- as.data.frame(colMeans(abadultAkey[sapply(abadultAkey, is.numeric)]))
names(new) <- "mean_expression"

newboth <- as.data.frame(colMeans(abadultAkeyPey[sapply(abadultAkeyPey, is.numeric)]))
names(newboth) <- "mean_expression"
newboth1 <- newboth %>% filter(newboth$mean_expression > 6.15)
nrow(newboth1)
p<-ggplot(newboth1, aes(x=rownames(newboth1), y=newboth1$mean_expression)) + 
  geom_dotplot(binaxis='y', stackdir='center')+theme(legend.position = "none")+labs(title="Genes in Akey and Pey",x="", y = "Mean expression")
p + coord_flip()
```

#Extracting gene expression data from Allen Brain Atlas - 5 stages
```{r}
##Loading dataset
data("dataset_5_stages")
#unique(dataset_5_stages$structure) #Checking structures

#(Perform enrichment on Genes from Deserts - Skip)
#input_hyper = data.frame(results$hgnc_symbol, is_candidate=1)
#res_devel = aba_enrich(input_hyper, dataset='5_stages')

#Selecting genes ID and structures present in dataset
id <- unique(dataset_5_stages$ensembl_gene_id)
st <- unique(dataset_5_stages$structure)
st_allen <- paste("Allen",st, sep=":") 
#Expression data for all structures and genes
ab <- get_expression(structure_ids=st_allen, gene_ids = id, dataset='5_stages')

#Converting data to a dataframe with useful format for later
list1 = vector(mode="list")
for (r in 1:length(ab)){
  ab[[r]] <- t(ab[[r]]) #transpose
  list1 <- get_name(colnames(ab[[r]])) #change Allen:XXXX to e.g. M1C_primary motor cortex, etc
  colnames(ab[[r]]) <- list1
  ab[[r]] <- as.data.frame(ab[[r]])
}
```

#SELECTION OF > q75 and generating file for genes within Akey
```{r}
#>q75 to select genes with high expression
q75 = vector(mode="list", length = length(ab))
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    q75[[i]][[h]] <- ab[[i]][h] %>% filter(ab[[i]][h] > quantile(ab[[i]][[h]], 0.75))
  }
}
#Changing row names ENSG to hgnc_symbol (via bioMart); ordering values based on expression
##ab is a list of 5 elements (i); each element or 'sublist' has 16 dataframes (h)
##q75[[1]][[3]][[1]] would denote  first sublist (1), third dataframe (3), column 1
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(q75[[i]][[h]]),mart=ensembl)
    q75[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(q75[[i]][[h]]), G_list$ensembl_gene_id)]
    q75[[i]][[h]] <-q75[[i]][[h]][order(q75[[i]][[h]][[1]], decreasing = TRUE), ]
  }
}
#Intersecting q75 (high expression) with Akey
akey075 = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    akey075[[i]][[h]] <- q75[[i]][[h]][q75[[i]][[h]][[2]] %in% results$hgnc_symbol,] #in Akey
    akey075[[i]][[h]] <-  akey075[[i]][[h]][!(is.na(akey075[[i]][[h]][[2]]) | akey075[[i]][[h]][[2]]==""), ] #Cleaning
  }
}
#Save results - Genes in deserts and pos. selection with high expression
for (i in 1:length(ab)){
  for (h in 1:length(akey075[[i]])){
    if (dim(akey075[[i]][[h]])[1] == 0) next #To skip empty dataframes
    write.xlsx(akey075[[i]][[h]], file="ABA_Akey_highExprq75.xlsx", sheetName=paste(toString(i),toString(names(akey075[[i]][[h]][1])), sep="_"), append = TRUE) #Add age number at the beginning of the sheet name
  }
}
```


#ABA Data - All genes that are present in Akey
```{r}
#Global data
ab1 <- ab
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
  write.xlsx(arrange(final_merge[i], desc(final_merge[i])), file="ABA_GenesAkey.xlsx", sheetName=names(final_merge[i]), append = TRUE)
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
ggsave(file="ABA_GenesAkey.pdf", pl1, width = 11.69, height = 8.27, units = "in")
```


#ABA Data - All genes that are present in both Akey and Pey
```{r}
#Global data
ab2 <- ab
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    ab2[[i]][[h]] <- as.data.frame(ab[[i]][h])
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(ab2[[i]][[h]]),mart=ensembl)
    ab2[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(ab2[[i]][[h]]), G_list$ensembl_gene_id)]
    ab2[[i]][[h]] <-ab2[[i]][[h]][order(ab2[[i]][[h]][[1]], decreasing = TRUE), ]
  }
}

aba_both = vector(mode="list", length = length(ab2))
for (i in 1:length(ab2)){
  for (h in 1:length(names(ab2[[i]]))){
    aba_both[[i]][[h]] <- ab2[[i]][[h]][ab2[[i]][[h]][[2]] %in% both$hgnc_symbol,] #in both
    aba_both[[i]][[h]] <-  aba_both[[i]][[h]][!(is.na(aba_both[[i]][[h]][[2]]) | aba_both[[i]][[h]][[2]]==""), ] #Cleaning
  }
}
#Reporting mean
mean2 = vector(mode="list", length = length(ab2))
for (i in 1:length(ab2)){
  for (h in 1:length(names(ab2[[i]]))){
    mean2[[i]][[h]] <- mean(aba_both[[i]][[h]][[1]])
    mean2[[i]][[h]] <- as.data.frame(mean2[[i]][[h]])
    names(mean2[[i]][[h]]) <- paste(names(aba_both[[i]][[h]][1]), sep='_')
  }
}

prenatal <- do.call("rbind", as.data.frame(mean2[[1]]))
colnames(prenatal) <- "prenatal"
infant <- do.call("rbind", as.data.frame(mean2[[2]]))
colnames(infant) <- "infant"
child <- do.call("rbind", as.data.frame(mean2[[3]]))
colnames(child) <- "child"
adolescent <- do.call("rbind", as.data.frame(mean2[[4]]))
colnames(adolescent) <- "adolescent"
adult <- do.call("rbind", as.data.frame(mean2[[5]]))
colnames(adult) <- "adult"

final_merge2 <- as.data.frame(cbind(prenatal, infant, child, adolescent, adult))
#Save results - need to be updated
library(xlsx)
for (i in 1:length(colnames(final_merge2))){
  write.xlsx(arrange(final_merge2[i], desc(final_merge2[i])), file="ABA_GenesAkeyPey.xlsx", sheetName=names(final_merge2[i]), append = TRUE)
}

#For plot:
final_merge2 <- tibble::rownames_to_column(final_merge2, "Structure")
final_merge2[[1]] <- sapply(strsplit(final_merge2[[1]], "_"), "[", 1)

a<-ggparcoord(final_merge2,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts and Pey")+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("mean  expression")
b <- ggparcoord(final_merge2,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("mean expression")
c<-ggparcoord(final_merge2,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("mean expression")

a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)
grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))

d<-ggparcoord(final_merge2,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Somato - Motor - Parietal - Aud Ctx")+scale_color_manual(values = c( "#00FF00", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#00FF00", "#ABABAB",  "#00FF00", "#ABABAB", "#ABABAB","#ABABAB", "#00FF00", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("mean expression")
d1<-arrangeGrob(d, left=textGrob("D"))

pl2 <- grid.arrange(a1, c1, d1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
ggsave(file="ABA_GenesAkeyPey.pdf", pl2, width = 11.69, height = 8.27, units = "in")

```

#ABA Data - All genes that are present in Pey
```{r}
#Global data
ab3 <- ab
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    ab3[[i]][[h]] <- as.data.frame(ab[[i]][h])
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(ab3[[i]][[h]]),mart=ensembl)
    ab3[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(ab3[[i]][[h]]), G_list$ensembl_gene_id)]
    ab3[[i]][[h]] <-ab3[[i]][[h]][order(ab3[[i]][[h]][[1]], decreasing = TRUE), ]
    
  }
}
aba_pey = vector(mode="list", length = length(ab3))
for (i in 1:length(ab3)){
  for (h in 1:length(names(ab3[[i]]))){
    aba_pey[[i]][[h]] <- ab3[[i]][[h]][ab3[[i]][[h]][[2]] %in% results_pey$hgnc_symbol,] #in Pey
    aba_pey[[i]][[h]] <-  aba_pey[[i]][[h]][!(is.na(aba_pey[[i]][[h]][[2]]) | aba_pey[[i]][[h]][[2]]==""), ] #Cleaning
  }
}
#Reporting mean
mean3 = vector(mode="list", length = length(ab3))
for (i in 1:length(ab3)){
  for (h in 1:length(names(ab3[[i]]))){
    mean3[[i]][[h]] <- mean(aba_pey[[i]][[h]][[1]])
    mean3[[i]][[h]] <- as.data.frame(mean3[[i]][[h]])
    names(mean3[[i]][[h]]) <- paste(names(aba_pey[[i]][[h]][1]), sep='_')
  }
}

prenatal <- do.call("rbind", as.data.frame(mean3[[1]]))
colnames(prenatal) <- "prenatal"
infant <- do.call("rbind", as.data.frame(mean3[[2]]))
colnames(infant) <- "infant"
child <- do.call("rbind", as.data.frame(mean3[[3]]))
colnames(child) <- "child"
adolescent <- do.call("rbind", as.data.frame(mean3[[4]]))
colnames(adolescent) <- "adolescent"
adult <- do.call("rbind", as.data.frame(mean3[[5]]))
colnames(adult) <- "adult"

final_merge3 <- as.data.frame(cbind(prenatal, infant, child, adolescent, adult))
#Save:
for (i in 1:length(colnames(final_merge3))){
  if (dim(final_merge3[i])[1] == 0) next
  write.xlsx(arrange(final_merge3[i], desc(final_merge3[i])), file="ABA_GenesPey.xlsx", sheetName=names(final_merge3[i]), append = TRUE)
}

#Plot
##For plot
final_merge3 <- tibble::rownames_to_column(final_merge3, "Structure")
final_merge3[[1]] <- sapply(strsplit(final_merge3[[1]], "_"), "[", 1)

a<-ggparcoord(final_merge3,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Pey")+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("mean  expression")
b <- ggparcoord(final_merge3,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("mean expression")
c<-ggparcoord(final_merge3,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("mean expression")
a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)

pl3 <- grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
ggsave(file="ABA_GenesPey.pdf", pl3, width = 11.69, height = 8.27, units = "in")
```

#ABA Data - All genes that are present in Rac
```{r}
#Global data
ab4 <- ab
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    ab4[[i]][[h]] <- as.data.frame(ab[[i]][h])
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(ab4[[i]][[h]]),mart=ensembl)
    ab4[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(ab4[[i]][[h]]), G_list$ensembl_gene_id)]
    ab4[[i]][[h]] <-ab4[[i]][[h]][order(ab4[[i]][[h]][[1]], decreasing = TRUE), ]
    
  }
}
aba_rac = vector(mode="list", length = length(ab4))
for (i in 1:length(ab4)){
  for (h in 1:length(names(ab4[[i]]))){
    aba_rac[[i]][[h]] <- ab4[[i]][[h]][ab4[[i]][[h]][[2]] %in% results_rac$hgnc_symbol,] #in Rac
    aba_rac[[i]][[h]] <-  aba_rac[[i]][[h]][!(is.na(aba_rac[[i]][[h]][[2]]) | aba_rac[[i]][[h]][[2]]==""), ] #Cleaning
  }
}
#Reporting mean
mean4 = vector(mode="list", length = length(ab4))
for (i in 1:length(ab4)){
  for (h in 1:length(names(ab4[[i]]))){
    mean4[[i]][[h]] <- mean(aba_rac[[i]][[h]][[1]])
    mean4[[i]][[h]] <- as.data.frame(mean4[[i]][[h]])
    names(mean4[[i]][[h]]) <- paste(names(aba_rac[[i]][[h]][1]), sep='_')
  }
}

prenatal <- do.call("rbind", as.data.frame(mean4[[1]]))
colnames(prenatal) <- "prenatal"
infant <- do.call("rbind", as.data.frame(mean4[[2]]))
colnames(infant) <- "infant"
child <- do.call("rbind", as.data.frame(mean4[[3]]))
colnames(child) <- "child"
adolescent <- do.call("rbind", as.data.frame(mean4[[4]]))
colnames(adolescent) <- "adolescent"
adult <- do.call("rbind", as.data.frame(mean4[[5]]))
colnames(adult) <- "adult"

final_merge4 <- as.data.frame(cbind(prenatal, infant, child, adolescent, adult))
#Save results_rac - need to be updated
for (i in 1:length(colnames(final_merge4))){
  if (dim(final_merge4[i])[1] == 0) next
  write.xlsx(arrange(final_merge4[i], desc(final_merge4[i])), file="ABA_GenesRac.xlsx", sheetName=names(final_merge4[i]), append = TRUE)
}

#Plot

##For plot
final_merge4 <- tibble::rownames_to_column(final_merge4, "Structure")
final_merge4[[1]] <- sapply(strsplit(final_merge4[[1]], "_"), "[", 1)

a<-ggparcoord(final_merge4,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Rac")+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("mean  expression")
b <- ggparcoord(final_merge4,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("mean expression")
c<-ggparcoord(final_merge4,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("mean expression")
a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))

grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)

pl4 <- grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
ggsave(file="ABA_GenesRac.pdf", pl4, width = 11.69, height = 8.27, units = "in")
```

#ABA Data - All genes that are present in both Akey and Rac
```{r}
#Global data
ab5 <- ab
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    ab5[[i]][[h]] <- as.data.frame(ab[[i]][h])
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(ab5[[i]][[h]]),mart=ensembl)
    ab5[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(ab5[[i]][[h]]), G_list$ensembl_gene_id)]
    ab5[[i]][[h]] <-ab5[[i]][[h]][order(ab5[[i]][[h]][[1]], decreasing = TRUE), ]
    
  }
}
aba_racAkey = vector(mode="list", length = length(ab5))
for (i in 1:length(ab5)){
  for (h in 1:length(names(ab5[[i]]))){
    aba_racAkey[[i]][[h]] <- ab5[[i]][[h]][ab5[[i]][[h]][[2]] %in% racAkey$hgnc_symbol,] #in RacAkey
    aba_racAkey[[i]][[h]] <-  aba_racAkey[[i]][[h]][!(is.na(aba_racAkey[[i]][[h]][[2]]) | aba_racAkey[[i]][[h]][[2]]==""), ] #Cleaning
  }
}
#Reporting mean
mean5 = vector(mode="list", length = length(ab5))
for (i in 1:length(ab5)){
  for (h in 1:length(names(ab5[[i]]))){
    mean5[[i]][[h]] <- mean(aba_racAkey[[i]][[h]][[1]])
    mean5[[i]][[h]] <- as.data.frame(mean5[[i]][[h]])
    names(mean5[[i]][[h]]) <- paste(names(aba_racAkey[[i]][[h]][1]), sep='_')
  }
}

prenatal <- do.call("rbind", as.data.frame(mean5[[1]]))
colnames(prenatal) <- "prenatal"
infant <- do.call("rbind", as.data.frame(mean5[[2]]))
colnames(infant) <- "infant"
child <- do.call("rbind", as.data.frame(mean5[[3]]))
colnames(child) <- "child"
adolescent <- do.call("rbind", as.data.frame(mean5[[4]]))
colnames(adolescent) <- "adolescent"
adult <- do.call("rbind", as.data.frame(mean5[[5]]))
colnames(adult) <- "adult"

final_merge5 <- as.data.frame(cbind(prenatal, infant, child, adolescent, adult))
#Save racAkey - need to be updated
for (i in 1:length(colnames(final_merge5))){
  if (dim(final_merge5[i])[1] == 0) next
  write.xlsx(arrange(final_merge5[i], desc(final_merge5[i])), file="ABA_GenesAkeyRac.xlsx", sheetName=names(final_merge5[i]), append = TRUE)
}

#Plot

##For plot
final_merge5 <- tibble::rownames_to_column(final_merge5, "Structure")
final_merge5[[1]] <- sapply(strsplit(final_merge5[[1]], "_"), "[", 1)

a<-ggparcoord(final_merge5,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Rac and Akey")+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("mean  expression")
b <- ggparcoord(final_merge5,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("mean expression")
c<-ggparcoord(final_merge5,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("mean expression")
a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)

pl5 <- grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
ggsave(file="ABA_GenesAkeyRac.pdf", pl5, width = 11.69, height = 8.27, units = "in")
```

#Table Cell Atlas - Science
```{sh}
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156793/suppl/GSE156793%5FS8%5FDE%5Fgene%5Fcells%2Ecsv%2Egz 
```

```{r}
library(dplyr)
df <- read.csv(file = "~/GSE156793_S8_DE_gene_cells.csv.gz")
df1 <- df[which(df$organ=='Cerebellum'), ]
write.csv(df1, file="CellAtlas_GSE156793_cerebellum.csv", row.names = FALSE)

df1$gene_short_name <- gsub("\\'", "", df1$gene_short_name)
df1subset <- df1[df1$gene_short_name %in% results$hgnc_symbol,]
df1subset <- df1subset[order(df1subset$max.cluster),]

df1subsetboth <- df1[df1$gene_short_name %in% both$hgnc_symbol,]
df1subsetboth <- df1subsetboth[order(df1subsetboth$max.cluster),]
df1subset <- df1subset[order(df1subset$max.cluster),]


write.csv(df1subset, file="CellAtlas_GSE156793_inAkey.csv", row.names = FALSE)
write.csv(df1subsetboth, file="CellAtlas_GSE156793_inAkeyPey.csv", row.names = FALSE)

#Whole dataset - AKEY
df$gene_short_name <- gsub("\\'", "", df$gene_short_name)
dfakey <- df[df$gene_short_name %in% results$hgnc_symbol,]
#dfakey <- dfakey[order(dfakey$max.cluster),]
#dfakey %>% group_by(organ) %>% mutate(mean1=mean(max.expr), mean2=mean(second.expr))
meanakey <- dfakey %>% group_by(organ) %>% summarize(Mean = mean(max.expr))
meanakey[order(meanakey$Mean, decreasing = TRUE),]
##Gene expression of genes in Akey - Preliminar test
pairwise.t.test(dfakey$max.expr, dfakey$organ, p.adjust.method = "BH")

#Whole dataset - AKEY PEY
df1subsetboth <- df[df$gene_short_name %in% both$hgnc_symbol,]
meanakeypey <- df1subsetboth %>% group_by(organ) %>% summarize(Mean = mean(max.expr), .groups = 'drop')
meanakeypey[order(meanakeypey$Mean, decreasing = TRUE),]

#Raw:
df  %>% group_by(organ) %>% summarize(Mean = mean(max.expr), .groups = 'drop')
```

#Friedman test - Genes in Akey:
```{r}
r <- data.frame(dfakey$organ, dfakey$gene_short_name, dfakey$max.expr)
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
library(tidyverse)
pp <-reduce(list_dfs, full_join, by="group") #Generating matrix for plot
pp[is.na(pp)] <- 0 #replacing null values
rownames(pp) <- pp[,1]
pp[,1] <- NULL
pp <- pp[c(2:15,1)] #Reordering columns for triangular matrix

#Plot
##Cell Atlas - Akey - Enrichment
library(Matrix)
library(reshape2)
library(ggplot2)
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
upper_tri <- get_upper_tri(as.matrix(pp))
melted <- melt(upper_tri, na.rm = TRUE)

pl_fr <- ggplot(data = melted, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#F0E442", high = "#0072B2", mid = "#0072B2", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
   name="p-value") +
  theme_classic()+xlab("")+ylab("")+ 
 theme(axis.text.x = element_text(angle = 45,hjust = 0))+scale_y_discrete(position = "right")+
 coord_fixed()+ coord_flip()
pl_fr
ggsave(file="CellAtlas_MeanExpr_heatmap.pdf", pl_fr,width = 11.69, height = 8.27, units = "in")

```


