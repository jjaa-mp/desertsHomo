#Packages
```{r}
#Packages
library(biomaRt)
library(ABAData)
library(ABAEnrichment)
library(data.table)
library(ggplot2)
library(dplyr)
library(GGally)
library(viridis)
library(ggpubr)
library(grid)
library(gridExtra)
library(lattice)
library(xlsx)
library(tidyr)
library(tidyverse)
library(DescTools)
library(reshape2)
library(dplyr)
library(Matrix)
library(pheatmap)
library(plyr)
library(rstatix)
library(matrixStats)
```

#Extracting gene names from coordinate regions via bioMart
```{r}
#Genes within Desert region coordinates
##Using hg19 genome
ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = 'hsapiens_gene_ensembl',
                        host = 'https://grch37.ensembl.org')
#Total number of protein coding genes in Ensembl
resP=getBM(attributes = c("hgnc_symbol","gene_biotype"),
              filters = "biotype",
              values = list(biotype="protein_coding"), mart = ensembl)
length(resP$hgnc_symbol[!(is.na(resP$hgnc_symbol) | resP$hgnc_symbol=="")])

##Akey Deserts coordinates (Ensembl 1-based):
filterlist <- c("1:105400000:120600000", "3:74100000:89300000", "7:106200000:123200000","8:49400000:66500000")
##Select only protein-coding from Akey
results=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=filterlist,biotype="protein_coding"), mart = ensembl)
results <- results[!duplicated(results$hgnc_symbol),]
#255
length(results$hgnc_symbol[!(is.na(results$hgnc_symbol) | results$hgnc_symbol=="")])
results <-  results[!(is.na(results$hgnc_symbol) | results$hgnc_symbol==""), ] #Cleaning


#ALTERNATIVE: 800 genes (protein coding genes plus other genes)
##results_alternative =getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
              #filters = c("chromosomal_region"),
              #values = list(chromosomal_region=filterlist), mart = ensembl)
##Save:
##write.csv(results, file="2020_Genes_in_Deserts.csv")

##Pey coordinates (Ensembl 1-based):
pey_coords <- read.delim("~/2020_pey_coords.bed", header=FALSE) #File with start 0-based
###Preparing bed file for input in bioMart
pey_coords$V1 <- gsub("chr", "\\1", pey_coords$V1)
pey_coords[2] <- pey_coords[2]+1 #Moving to 1-based for Ensembl
df <- paste(pey_coords$V1, pey_coords$V2, pey_coords$V3, sep = ":")
results_pey=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=df,biotype="protein_coding"), mart = ensembl)
results_pey <- results_pey[!duplicated(results_pey$hgnc_symbol),]
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
results_rac <- results_rac[!duplicated(results_rac$hgnc_symbol),]
###Genes within Akey and Rac - RESULT
racAkey <- results_rac[results_rac$hgnc_symbol %in% results$hgnc_symbol,]
racAkey <- racAkey[!(is.na(racAkey$hgnc_symbol) | racAkey$hgnc_symbol==""), ]
```

```{r}
library(ABAEnrichment)

#AKEY
input_hyper0 = data.frame(results$hgnc_symbol, is_candidate=1)
resABA0<-aba_enrich(input_hyper0, dataset="5_stages")
resABA00<-aba_enrich(input_hyper0, dataset="dev_effect")

#PEY
input_hyper1 = data.frame(both$hgnc_symbol, is_candidate=1)
resABA1<-aba_enrich(input_hyper1, dataset="5_stages")
resABA11<-aba_enrich(input_hyper1, dataset="dev_effect")


#Replicating Vernot
##desert from both nean and deni
listvernot <- c("1:104000000:114900000", "3:76500000:90500000", "7:113600000:124700000","8:54500000:65400000")
##Select only protein-coding from Akey
vernot=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=listvernot,biotype="protein_coding"), mart = ensembl)
vernot <- vernot[!duplicated(vernot$hgnc_symbol),]
vernot <-  vernot[!(is.na(vernot$hgnc_symbol) | vernot$hgnc_symbol==""), ]

input_vernot = data.frame(vernot$hgnc_symbol, is_candidate=1)
resvernot<-aba_enrich(input_vernot, dataset="5_stages")
resvernot2<-aba_enrich(input_vernot, dataset="dev_effect")


listvernot2 <- c("1:102200000:114900000","2:201100000:211500000", "3:76500000:90500000", "7:106300000:124700000","8:53900000:66000000", "18:25000000:41800000")
##Select only protein-coding from Akey
vernotgenes2=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=listvernot2,biotype="protein_coding"), mart = ensembl)
vernotres2 <- vernotgenes2[!duplicated(vernotgenes2$hgnc_symbol),]

input_vernot2 = data.frame(vernotres2$hgnc_symbol, is_candidate=1)
vernotNean<-aba_enrich(input_vernot2, dataset="5_stages")
vernotNean2<-aba_enrich(input_vernot2, dataset="dev_effect")

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

#ABA Data - All genes that are present in Akey - LOG NORMALIZE
```{r}
ab1 = vector(mode="list", length = length(ab))
 for (i in 1:length(ab)){
   for (h in 1:length(names(ab[[i]]))){
     ab1[[i]][[h]] <- ab[[i]][h]
   }
 }
 
for (i in 1:length(ab)){
   for (h in 1:length(names(ab[[i]]))){
     G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(ab1[[i]][[h]]),mart=ensembl)
     ab1[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(ab1[[i]][[h]]), G_list$ensembl_gene_id)]
     ab1[[i]][[h]] <-ab1[[i]][[h]][order(ab1[[i]][[h]][[1]], decreasing = TRUE), ]
   }
}

#Intersecting q10 (high expression) with Akey
akey10 = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    akey10[[i]][[h]] <- ab1[[i]][[h]][ab1[[i]][[h]][[2]] %in% results$hgnc_symbol,] #in Akey
    akey10[[i]][[h]] <-  akey10[[i]][[h]][!(is.na(akey10[[i]][[h]][[2]]) | akey10[[i]][[h]][[2]]==""), ]
    akey10[[i]][[h]][1] <- log2(akey10[[i]][[h]][1]+1)
    #Cleaning
  }
}

#Reporting mean
mean1 = vector(mode="list", length = length(ab))
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    mean1[[i]][[h]] <- median(akey10[[i]][[h]][[1]]) #mean
    mean1[[i]][[h]] <- as.data.frame(mean1[[i]][[h]])
    names(mean1[[i]][[h]]) <- paste(names(akey10[[i]][[h]][1]), sep='_')
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
#write.xlsx(final_merge,row.names = TRUE, file = "~/raul_tesina/1.data/ABAData_AkeyPeyRac_log2/ABA_GenesAkey_log2.xlsx")

#Plot
final_merge <- tibble::rownames_to_column(final_merge, "Structure")
final_merge[[1]] <- sapply(strsplit(final_merge[[1]], "_"), "[", 1)

a<-ggparcoord(final_merge,
              columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts")+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("expression")
b <- ggparcoord(final_merge,
                columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
c<-ggparcoord(final_merge,
              columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("expression")
a1<-arrangeGrob(a, left=textGrob("A")) 
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)

pl1abak <-grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
##ggsave(file="~/raul_tesina/2.plots/ABAData_AkeyPey_log2_median/ABA_GenesAkey_log_median.pdf", pl1abak, width = 11.69, height = 8.27, units = "in")

#STAT
dfa <- lapply(akey10, melt)
for (i in 1:length(dfa)){ 
   dfa[[i]][4] <- i #in Akey
 }
dfa1 <- ldply(dfa, data.frame)
colnames(dfa1)[colnames(dfa1) == "L1"] <- c("stage")
dfa1$stage <- as.factor(dfa1$stage)
dfa1 <- as_tibble(dfa1)
dfa1 <- dfa1 %>% mutate(variable=as.character(variable))
# 
# dfa_s1 <- subset(dfa1, dfa1$stage==1) %>% droplevels()
# dfa_s2 <- subset(dfa1, dfa1$stage==2) %>% droplevels()
# dfa_s3 <- subset(dfa1, dfa1$stage==3) %>% droplevels()
# dfa_s4 <- subset(dfa1, dfa1$stage==4) %>% droplevels()
# dfa_s5 <- subset(dfa1, dfa1$stage==5) %>% droplevels()
# dfa_s1 <- dfa_s1[dfa_s1$gene_name %in% names(which(table(dfa_s1$gene_name) > 15)), ]
# dfa_s2 <- dfa_s2[dfa_s2$gene_name %in% names(which(table(dfa_s2$gene_name) > 15)), ]
# dfa_s3 <- dfa_s3[dfa_s3$gene_name %in% names(which(table(dfa_s3$gene_name) > 15)), ]
# dfa_s4 <- dfa_s4[dfa_s4$gene_name %in% names(which(table(dfa_s4$gene_name) > 15)), ]
# dfa_s5 <- dfa_s5[dfa_s5$gene_name %in% names(which(table(dfa_s5$gene_name) > 15)), ]
 
finaldf <- dfa1[dfa1$gene_name %in% names(which(table(dfa1$gene_name) > 79)), ]
#ggqqplot(finaldf, "value", facet.by = "stage")
#qq1 <- ggqqplot(finaldf, "value", facet.by = "stage")
###ggsave(file="~/raul_tesina/2.plots/ABAData_AkeyPeyRac_log2/ABA_GenesAkey_qqplot.pdf", qq1, width = 11.69, height = 8.27, units = "in")


oneway1 <- finaldf %>%
  group_by(stage) %>%
  anova_test(dv = value, wid = gene_name, within = variable) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
oneway1

#pairwise comparisons
pw1 <- finaldf %>%
  group_by(stage) %>%
  pairwise_t_test(
    value ~ variable, paired = TRUE,
    p.adjust.method = "BH"
    )
struc_pw <- pw1 %>% filter(p.adj.signif != "ns")
##write.csv(struc_pw, file="~/raul_tesina/1.data/ABAData_AkeyPeyRac_log2/ABA_GenesAkey_log2_anova_Structures_pairwise_significant.csv", row.names = FALSE)

oneway2 <- finaldf %>%
  group_by(variable) %>%
  anova_test(dv = value, wid = gene_name, within = stage) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "BH")
oneway2 #one is significant

#pairwise comparisons
pw2 <- finaldf %>%
  group_by(variable) %>%
  pairwise_t_test(
    value ~ stage, paired = TRUE,
    p.adjust.method = "BH"
    )
stages <- pw2 %>% filter(p.adj.signif != "ns")
##write.csv(stages, file="~/raul_tesina/1.data/ABAData_AkeyPeyRac_log2/ABA_GenesAkey_log2_anova__Stages_significant.csv", row.names = FALSE)

```

#ABA Data - All genes that are present in both Akey and Pey LOG NORMALIZE
```{r}
#Global data
ab2 = vector(mode="list", length = length(ab))
 for (i in 1:length(ab)){
   for (h in 1:length(names(ab[[i]]))){
     ab2[[i]][[h]] <- ab[[i]][h]
   }
 }
 
for (i in 1:length(ab)){
   for (h in 1:length(names(ab[[i]]))){
     G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(ab2[[i]][[h]]),mart=ensembl)
     ab2[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(ab2[[i]][[h]]), G_list$ensembl_gene_id)]
     ab2[[i]][[h]] <-ab2[[i]][[h]][order(ab2[[i]][[h]][[1]], decreasing = TRUE), ]
   }
}
akeypey10 = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    akeypey10[[i]][[h]] <- ab2[[i]][[h]][ab2[[i]][[h]][[2]] %in% both$hgnc_symbol,] #in Akey
    akeypey10[[i]][[h]] <-  akeypey10[[i]][[h]][!(is.na(akeypey10[[i]][[h]][[2]]) | akeypey10[[i]][[h]][[2]]==""), ]
    akeypey10[[i]][[h]][1] <- log2(akeypey10[[i]][[h]][1]+1)
#Cleaning
  }
}

mean2 = vector(mode="list", length = length(ab))
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    mean2[[i]][[h]] <- median(akeypey10[[i]][[h]][[1]]) #mean
    mean2[[i]][[h]] <- as.data.frame(mean2[[i]][[h]])
    names(mean2[[i]][[h]]) <- paste(names(akeypey10[[i]][[h]][1]), sep='_')
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

final_mergelog <- as.data.frame(cbind(prenatal, infant, child, adolescent, adult))
#write.xlsx(final_mergelog,row.names = TRUE, file = "~/raul_tesina/1.data/ABAData_AkeyPeyRac_log2/ABA_GenesAkeyPey_log2.xlsx")

#For plot:
final_mergelog <- tibble::rownames_to_column(final_mergelog, "Structure")
final_mergelog[[1]] <- sapply(strsplit(final_mergelog[[1]], "_"), "[", 1)

a<-ggparcoord(final_mergelog,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts and Pey")+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("expression")
b <- ggparcoord(final_mergelog,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum (red) & MD (green)")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#2ca25f", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
c<-ggparcoord(final_mergelog,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum (blue) and HIP (orange)")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#fdae6b", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("expression")

a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)
grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))

d<-ggparcoord(final_mergelog,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Somato - Motor - Parietal - Aud Ctx")+scale_color_manual(values = c( "#00FF00", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#00FF00", "#ABABAB",  "#00FF00", "#ABABAB", "#ABABAB","#ABABAB", "#00FF00", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("expression")
d1<-arrangeGrob(d, left=textGrob("D"))

plakap <- grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
##ggsave(file="~/raul_tesina/2.plots/ABAData_AkeyPey_log2_median/ABA_GenesAkeyPey_log_median.pdf", plakap, width = 11.69, height = 8.27, units = "in")

# STAT
dfa <- lapply(akeypey10, melt)
for (i in 1:length(dfa)){ 
   dfa[[i]][4] <- i #in Akey
 }
dfa2 <- ldply(dfa, data.frame)
colnames(dfa2)[colnames(dfa2) == "L1"] <- c("stage")
dfa2$stage <- as.factor(dfa2$stage)
dfa2 <- as_tibble(dfa2)
dfa2 <- dfa2 %>% mutate(variable=as.character(variable))

# dfa_s1 <- subset(dfa2, dfa2$stage==1) %>% droplevels()
# dfa_s2 <- subset(dfa2, dfa2$stage==2) %>% droplevels()
# dfa_s3 <- subset(dfa2, dfa2$stage==3) %>% droplevels()
# dfa_s4 <- subset(dfa2, dfa2$stage==4) %>% droplevels()
# dfa_s5 <- subset(dfa2, dfa2$stage==5) %>% droplevels()
# 
# dfa_s1 <- dfa_s1[dfa_s1$gene_name %in% names(which(table(dfa_s1$gene_name) > 15)), ]
# dfa_s2 <- dfa_s2[dfa_s2$gene_name %in% names(which(table(dfa_s2$gene_name) > 15)), ]
# dfa_s3 <- dfa_s3[dfa_s3$gene_name %in% names(which(table(dfa_s3$gene_name) > 15)), ]
# dfa_s4 <- dfa_s4[dfa_s4$gene_name %in% names(which(table(dfa_s4$gene_name) > 15)), ]
# dfa_s5 <- dfa_s5[dfa_s5$gene_name %in% names(which(table(dfa_s5$gene_name) > 15)), ]


finaldf2 <- dfa2[dfa2$gene_name %in% names(which(table(dfa2$gene_name) > 79)), ]
```

#ABA Data - All genes (raw data) LOG NORMALIZE
```{r}
#ab
aba_all = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    aba_all[[i]][[h]] <- q10[[i]][[h]][q10[[i]][[h]][[2]] %in% resP$hgnc_symbol,]
    aba_all[[i]][[h]] <-  aba_all[[i]][[h]][!(is.na(aba_all[[i]][[h]][[2]]) | aba_all[[i]][[h]][[2]]==""), ]
    aba_all[[i]][[h]][1] <- log2(aba_all[[i]][[h]][1])
#Cleaning
  }
}

mean0 = vector(mode="list", length = length(ab))
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    mean0[[i]][[h]] <- mean(aba_all[[i]][[h]][[1]])
    mean0[[i]][[h]] <- as.data.frame(mean0[[i]][[h]])
    names(mean0[[i]][[h]]) <- paste(names(aba_all[[i]][[h]][1]), sep='_')
  }
}

prenatal <- do.call("rbind", as.data.frame(mean0[[1]]))
colnames(prenatal) <- "prenatal"
infant <- do.call("rbind", as.data.frame(mean0[[2]]))
colnames(infant) <- "infant"
child <- do.call("rbind", as.data.frame(mean0[[3]]))
colnames(child) <- "child"
adolescent <- do.call("rbind", as.data.frame(mean0[[4]]))
colnames(adolescent) <- "adolescent"
adult <- do.call("rbind", as.data.frame(mean0[[5]]))
colnames(adult) <- "adult"

final_merge0 <- as.data.frame(cbind(prenatal, infant, child, adolescent, adult))
#For plot
final_merge0 <- tibble::rownames_to_column(final_merge0, "Structure")
final_merge0[[1]] <- sapply(strsplit(final_merge0[[1]], "_"), "[", 1)

a<-ggparcoord(final_merge0,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes Raw Data")+scale_color_viridis(discrete=TRUE)+theme(plot.title = element_text(size=10))+xlab("")+ylab("mean  expression")
b <- ggparcoord(final_merge0,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
c<-ggparcoord(final_merge0,
    columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("expression")
a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)


pl3 <-grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
###ggsave(file="~/raul_tesina/2.plots/ABAData_AkeyPeyRac_log2/ABA_GenesAll_log.pdf", pl3, width = 11.69, height = 8.27, units = "in")


dfa <- lapply(aba_all, melt)
for (i in 1:length(dfa)){ 
   dfa[[i]][4] <- i #in Akey
 }
dfa1 <- ldply(dfa, data.frame)
colnames(dfa1)[colnames(dfa1) == "L1"] <- c("stage")
dfa1$stage <- as.factor(dfa1$stage)
dfa1 <- as_tibble(dfa1)
dfa1 <- dfa1 %>% mutate(variable=as.character(variable))
finaldf0 <- dfa1[dfa1$gene_name %in% names(which(table(dfa1$gene_name) > 79)), ]


```

#ABA - STATS
```{R}
#Friedman test per each stage independently.
##Post hoc via Conover
###AKEY & PEY
finaldf2
aba_akeypey_st <- lapply(split(finaldf2, finaldf2$stage), function(x) {friedman_test(value ~ variable | gene_name, data=x)}) #All except stage 1

#Friedman test per each stage independently.
###ABA - AKEY alone
finaldf
aba_akey_st <- lapply(split(finaldf, finaldf$stage), function(x) {friedman_test(value ~ variable | gene_name, data=x)}) #All significant
```

#Table Cell Atlas - Science
```{sh}
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156793/suppl/GSE156793%5FS8%5FDE%5Fgene%5Fcells%2Ecsv%2Egz 
```


```{r}
df <- read.csv(file = "~/GSE156793_S8_DE_gene_cells.csv.gz")

#Whole dataset  in AKEY
df <- df %>% mutate(max.expr = log2(max.expr))
df$gene_short_name <- gsub("\\'", "", df$gene_short_name)
dfakey <- df %>% filter(gene_short_name %in% results$hgnc_symbol)
meanakey <- dfakey %>% group_by(organ) %>% dplyr::summarize(Mean = mean(max.expr))
meanakey[order(meanakey$Mean, decreasing = TRUE),]





#filtering for paired test
meanakeytest <- dfakey[dfakey$gene_short_name %in% names(which(table(dfakey$gene_short_name) > 14)), ]
meanakeytest <- meanakeytest %>% arrange(meanakeytest$gene_short_name)
pp <- pairwise.t.test(meanakeytest$max.expr, meanakeytest$organ, data=meanakeytest,  p.adjust.method = "BH", paired = TRUE)

get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
lower_tri <- get_lower_tri(as.matrix(pp$p.value))
melted <- melt(lower_tri, na.rm = TRUE)

pl_fr <- ggplot(data = melted, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#F0E442", high = "#0072B2", mid = "#CC79A7", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
   name="p-value") +
  theme_classic()+xlab("")+ylab("")+ 
 theme(axis.text.x = element_text(angle = 45,hjust = 0))+scale_x_discrete(position = "top")+
 coord_fixed()
pl_fr
##ggsave(file="~/raul_tesina/2.plots/cellatlas_meanexpr_log_heatmap/CellAtlas_MeanExprAkey_log_heatmap.pdf", pl_fr,width = 11.69, height = 8.27, units = "in")


#Whole dataset in AKEY & PEY
dfsubsetboth <- df %>% filter(gene_short_name %in% both$hgnc_symbol)
meanakeypey <- dfsubsetboth %>% group_by(organ) %>% dplyr::summarize(Mean = mean(max.expr))
meanakeypey[order(meanakeypey$Mean, decreasing = TRUE),]

#filtering for paired test
meanakeypeytest <- dfsubsetboth[dfsubsetboth$gene_short_name %in% names(which(table(dfsubsetboth$gene_short_name) > 14)), ]
meanakeypeytest <- meanakeypeytest %>% arrange(meanakeypeytest$gene_short_name)

pp1<-pairwise.t.test(meanakeypeytest$max.expr, meanakeypeytest$organ, p.adjust.method = "BH", paired = TRUE)

lower_tri <- get_lower_tri(as.matrix(pp1$p.value))
melted1 <- melt(lower_tri, na.rm = TRUE)

pl_fr1 <- ggplot(data = melted1, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#F0E442", high = "#0072B2", mid = "#CC79A7", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
   name="p-value") +
  theme_classic()+xlab("")+ylab("")+ 
 theme(axis.text.x = element_text(angle = 45,hjust = 0))+scale_x_discrete(position = "top")+
 coord_fixed()
pl_fr1
##ggsave(file="~/raul_tesina/2.plots/cellatlas_meanexpr_log_heatmap/CellAtlas_MeanExprAkeyPey_log_heatmap.pdf", pl_fr1,width = 11.69, height = 8.27, units = "in")


#Raw:
rawmean <- df  %>% group_by(organ) %>% dplyr::summarize(Mean = mean(max.expr))
rawmean[order(rawmean$Mean, decreasing = TRUE),]

df1 <- df[which(df$organ=='Cerebellum'), ]
brain <- df[which(df$organ=='Cerebrum'), ]
##write.csv(df1, file="CellAtlas_GSE156793_cerebellum.csv", row.names = FALSE)

# CBL - AKEY
df1$gene_short_name <- gsub("\\'", "", df1$gene_short_name)
df1subset <- df1[df1$gene_short_name %in% results$hgnc_symbol,]
df1subset <- df1subset[order(df1subset$max.cluster),]
# CBL - AKEY & PEY
df1subsetboth <- df1[df1$gene_short_name %in% both$hgnc_symbol,]
df1subsetboth <- df1subsetboth[order(df1subsetboth$max.cluster),]
df1subset <- df1subset[order(df1subset$max.cluster),]
#Cerebellum FILES:
##write.csv(df1subset, file="CellAtlas_GSE156793_CBL_inAkey.csv", row.names = FALSE)
##write.csv(df1subsetboth, file="CellAtlas_GSE156793_CBL_inAkeyPey.csv", row.names = FALSE)
```

#Friedman test - Genes in Akey:
```{r}
#Checking normality

r <- data.frame(dfakey$organ, dfakey$gene_short_name, dfakey$max.expr)
names(r) <- c("x", "y", "z")
rr <- pivot_wider(r, names_from = x, values_from = z)
rr <- rr[complete.cases(rr),]
rr3 <- as.matrix(rr)
rr3 <-rr3[,-1]
friedman.test(rr3)

list1 <- vector(mode = "list")
for (i in 1:length(colnames(rr3))){
  list1[[i]] <- as.numeric(c(rr3[,i]))
}
res_F <- DunnettTest(list1, control = c(1:15))
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

#Plot
##Cell Atlas - Akey - Enrichment

  # Get upper triangle of the correlation matrix
get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
}
lower_tri <- get_lower_tri(as.matrix(pp$p.value))
melted <- melt(lower_tri, na.rm = TRUE)

pl_fr <- ggplot(data = melted, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#F0E442", high = "#0072B2", mid = "#CC79A7", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
   name="p-value") +
  theme_classic()+xlab("")+ylab("")+ 
 theme(axis.text.x = element_text(angle = 45,hjust = 0))+scale_y_discrete(position = "right")+
 coord_fixed()+ coord_flip()
pl_fr
###ggsave(file="CellAtlas_MeanExprAkey_heatmap.pdf", pl_fr,width = 11.69, height = 8.27, units = "in")

```

#Friedman test - Genes in both Akey and Pey:
```{r}
ggplot(s, aes(z)) + geom_histogram()

ggplot(s[which(s$z > 0 & s$z < 500),], aes(z)) +
    geom_histogram()

s <- data.frame(dfsubsetboth$organ, dfsubsetboth$gene_short_name, dfsubsetboth$max.expr)
names(s) <- c("x", "y", "z")
ss <- pivot_wider(s, names_from = x, values_from = z)
ss <- ss[complete.cases(ss),]
ss3 <- as.matrix(ss)
ss3 <-ss3[,-1]
friedman.test(ss3)

list2 <- vector(mode = "list")
for (i in 1:length(colnames(ss3))){
  list2[[i]] <- as.numeric(c(ss3[,i]))
}
res_F2 <- DunnettTest(list2, control = c(1:15)) #CAN ALSO TRY SNK Test
organs <- colnames(ss3)

list_dfs2 <- vector(mode="list")
for (i in 1:length(res_F2)){
  colnames(res_F2[[i]]) <- c("diff", "lwr.ci", "upr.ci", colnames(ss3)[i]) #control group: colnames(ss3)[i]
  rownames(res_F2[[i]]) <- organs[-i] #comparison groups in rownames: all except i
  list_dfs2[[i]] <- as.data.frame(res_F2[[i]])
  list_dfs2[[i]] <- list_dfs2[[i]][4]
  list_dfs2[[i]] <-  cbind(list_dfs2[[i]], group=rownames(list_dfs2[[i]]))
  list_dfs2[[i]] <- list_dfs2[[i]][c(2,1)] #ordering columns for merging on index 'group'
}

pp2 <-reduce(list_dfs2, full_join, by="group") #Generating matrix for plot
pp2[is.na(pp2)] <- 2 #replacing null values
rownames(pp2) <- pp2[,1]
pp2[,1] <- NULL
pp2 <- pp2[c(2:15,1)] #Reordering columns for triangular matrix

#Plot
##Cell Atlas - Akey - Enrichment

# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
cormat[lower.tri(cormat)]<- NA
return(cormat)
}
upper_tri2 <- get_upper_tri(as.matrix(pp2))
melted2 <- melt(upper_tri2, na.rm = TRUE)

pl_fr2 <- ggplot(data = melted2, aes(Var2, Var1, fill = value))+
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "#F0E442", high = "#0072B2", mid = "#CC79A7", 
   midpoint = 0.5, limit = c(0,1), space = "Lab", 
   name="p-value") +
  theme_classic()+xlab("")+ylab("")+ 
 theme(axis.text.x = element_text(angle = 45,hjust = 0))+scale_y_discrete(position = "right")+
 coord_fixed()+ coord_flip()
pl_fr2
###ggsave(file="CellAtlas_MeanExprAkeyPey_heatmap.pdf", pl_fr2,width = 11.69, height = 8.27, units = "in")
```

#CELL ATLAS - CELL TYPE ENRICHMENT - AKEY
```{r}
#CEREBRUM - Hyper geom
brain$gene_short_name <- gsub("\\'", "", brain$gene_short_name)
brain_akey <- brain[brain$gene_short_name %in% results$hgnc_symbol,]
ctb <- data.frame(brain_akey$max.cluster, brain_akey$gene_short_name, brain_akey$max.expr)
names(ctb) <- c("x", "y", "z")
ctb1 <- pivot_wider(ctb, names_from = x, values_from = z)
colSums(!is.na(ctb1))#sample size 228
z0 <- as.data.frame(colSums(!is.na(ctb1[,-1]))) #removing y
names(z0) <- "Akey genes"
setDT(z0, keep.rownames = "Var1")
z0 <- as.data.frame(z0)
#testing
gr <- data.frame(brain$max.cluster, brain$gene_short_name, brain$max.expr)
names(gr) <- c("x", "y", "z")
gr1 <- as.data.frame(table(gr$x))
gr1<-gr1[which(gr1$Freq != 0),]
names(gr1) <- c("Var1", "Total Cerebrum DEG")
finaltable <- merge(gr1, z0, by="Var1")
##Hyper geom:
finaltable$pvalue <- phyper(finaltable$`Akey genes`-1, 255, 19238-255, finaltable$`Total Cerebrum DEG`, lower.tail = FALSE, log.p = FALSE)

#CBL - Hyper geom
df1$gene_short_name <- gsub("\\'",  "", df1$gene_short_name)
cbl_akey <- df1[df1$gene_short_name %in% results$hgnc_symbol,]
cblt <- data.frame(cbl_akey$max.cluster, cbl_akey$gene_short_name, cbl_akey$max.expr)
names(cblt) <- c("x", "y", "z")
cblt1 <- pivot_wider(cblt, names_from = x, values_from = z)
colSums(!is.na(cblt1))
cblt1
colSums(!is.na(cblt1))
z1 <- as.data.frame(colSums(!is.na(cblt1[,-1]))) #removing y
names(z1) <- "Akey genes"
setDT(z1, keep.rownames = "Var1")
z1 <- as.data.frame(z1)
#testing
fi1 <- data.frame(df1$max.cluster, df1$gene_short_name, df1$max.expr)
names(fi1) <- c("x", "y", "z")
fi1 <- as.data.frame(table(fi1$x))
fi1<-fi1[which(fi1$Freq != 0),]
#testing
fi1 <- data.frame(df1$max.cluster, df1$gene_short_name, df1$max.expr)
names(fi1) <- c("x", "y", "z")
fi1 <- as.data.frame(table(fi1$x))
fi1<-fi1[which(fi1$Freq != 0),]
names(fi1) <- c("Var1", "Total CBL DEG")
finaltable2 <- merge(fi1, z1, by="Var1")
finaltable2
finaltable2$pvalue <- phyper(finaltable2$`Akey genes`-1, 255, 19238-255, finaltable2$`Total CBL DEG`, lower.tail = FALSE, log.p = FALSE)


#CEREBRUM - AKEY & PEY
#brain$gene_short_name <- gsub("\\'", "", brain$gene_short_name)
brain_akeypey <- brain[brain$gene_short_name %in% both$hgnc_symbol,]
brain_akeypey1 <- data.frame(brain_akeypey$max.cluster, brain_akeypey$gene_short_name, brain_akeypey$max.expr)
names(brain_akeypey1) <- c("x", "y", "z")
brain_akeypey2 <- pivot_wider(brain_akeypey1, names_from = x, values_from = z)
colSums(!is.na(brain_akeypey2))
#CEREBELLUM - AKEY & PEY
#df1$gene_short_name <- gsub("\\'", "", df1$gene_short_name)
cbl_akeypey <- df1[df1$gene_short_name %in% both$hgnc_symbol,]
cbl_akeypey1 <- data.frame(cbl_akeypey$max.cluster, cbl_akeypey$gene_short_name, cbl_akeypey$max.expr)
names(cbl_akeypey1) <- c("x", "y", "z")
cbl_akeypey2 <- pivot_wider(cbl_akeypey1, names_from = x, values_from = z)
colSums(!is.na(cbl_akeypey2))
```

#CELL ATLAS - CELL TYPE ENRICHMENT - AKEYPEY (low numbers; no significant)
```{r}
# #CEREBRUM - Hyper geom
# brain$gene_short_name <- gsub("\\'", "", brain$gene_short_name)
# brain_akeypey <- brain[brain$gene_short_name %in% both$hgnc_symbol,]
# ctb <- data.frame(brain_akeypey$max.cluster, brain_akeypey$gene_short_name, brain_akeypey$max.expr)
# names(ctb) <- c("x", "y", "z")
# ctb1 <- pivot_wider(ctb, names_from = x, values_from = z)
# colSums(!is.na(ctb1))#sample size 228
# z0 <- as.data.frame(colSums(!is.na(ctb1[,-1]))) #removing y
# names(z0) <- "akeypey genes"
# setDT(z0, keep.rownames = "Var1")
# z0 <- as.data.frame(z0)
# #testing
# gr <- data.frame(brain$max.cluster, brain$gene_short_name, brain$max.expr)
# names(gr) <- c("x", "y", "z")
# gr1 <- as.data.frame(table(gr$x))
# gr1<-gr1[which(gr1$Freq != 0),]
# names(gr1) <- c("Var1", "Total Cerebrum DEG")
# finaltable <- merge(gr1, z0, by="Var1")
# ##Hyper geom:
# finaltable$pvalue <- phyper(finaltable$`akeypey genes`-1, 255, 19238-255, finaltable$`Total Cerebrum DEG`, lower.tail = FALSE, log.p = FALSE)
# 
# #CBL - Hyper geom
# df1$gene_short_name <- gsub("\\'",  "", df1$gene_short_name)
# cbl_akeypey <- df1[df1$gene_short_name %in% both$hgnc_symbol,]
# cblt <- data.frame(cbl_akeypey$max.cluster, cbl_akeypey$gene_short_name, cbl_akeypey$max.expr)
# names(cblt) <- c("x", "y", "z")
# cblt1 <- pivot_wider(cblt, names_from = x, values_from = z)
# colSums(!is.na(cblt1))
# cblt1
# colSums(!is.na(cblt1))
# z1 <- as.data.frame(colSums(!is.na(cblt1[,-1]))) #removing y
# names(z1) <- "akeypey genes"
# setDT(z1, keep.rownames = "Var1")
# z1 <- as.data.frame(z1)
# #testing
# fi1 <- data.frame(df1$max.cluster, df1$gene_short_name, df1$max.expr)
# names(fi1) <- c("x", "y", "z")
# fi1 <- as.data.frame(table(fi1$x))
# fi1<-fi1[which(fi1$Freq != 0),]
# #testing
# fi1 <- data.frame(df1$max.cluster, df1$gene_short_name, df1$max.expr)
# names(fi1) <- c("x", "y", "z")
# fi1 <- as.data.frame(table(fi1$x))
# fi1<-fi1[which(fi1$Freq != 0),]
# names(fi1) <- c("Var1", "Total CBL DEG")
# finaltable2 <- merge(fi1, z1, by="Var1")
# finaltable2
# finaltable2$pvalue <- phyper(finaltable2$`akeypey genes`-1, 255, 19238-255, finaltable2$`Total CBL DEG`, lower.tail = FALSE, log.p = FALSE)
# 
# 
# #EYE - Hyper geom
# eye <- df[which(df$organ=='Eye'), ]
# eye$gene_short_name <- gsub("\\'",  "", eye$gene_short_name)
# eye_akeypey <- eye[eye$gene_short_name %in% both$hgnc_symbol,]
# eyelt <- data.frame(eye_akeypey$max.cluster, eye_akeypey$gene_short_name, eye_akeypey$max.expr)
# names(eyelt) <- c("x", "y", "z")
# eyelt <- pivot_wider(eyelt, names_from = x, values_from = z)
# colSums(!is.na(eyelt)) #Sample size 215
# z2 <- as.data.frame(colSums(!is.na(eyelt[,-1]))) #removing y
# names(z2) <- "akeypey genes"
# setDT(z2, keep.rownames = "Var1")
# z2 <- as.data.frame(z2)
# #testing
# fi2 <- data.frame(eye$max.cluster, eye$gene_short_name, eye$max.expr)
# names(fi2) <- c("x", "y", "z")
# fi2 <- as.data.frame(table(fi2$x))
# fi2<-fi2[which(fi2$Freq != 0),]
# names(fi2) <- c("Var1", "Total EYE DEG")
# finaltable3 <- merge(fi2, z2, by="Var1")
# #hyper geom test
finaltable3$pvalue <- phyper(finaltable3$`akeypey genes`-1, 255, 19238-255, finaltable3$`Total EYE DEG`, lower.tail = FALSE, log.p = FALSE)
```

#Raw Sestan
```{r}
mRNAseqData=read.table("~/tmp_psychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
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
metadatamRNAseq=read_xlsx("~/tmp_psychENCODE/mRNA-seq_QC.xlsx",skip = 3)
modMetadatamRNAseq=na.omit(metadatamRNAseq)
modMetadatamRNAseq=modMetadatamRNAseq %>% select(1:3)
modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)
finalrawSestan=merge(modMetadatamRNAseq,lograwsestan1,by=c("Braincode", "Regioncode"))
finalraw2 <- finalrawSestan
rownames(finalraw2)<- do.call(paste,c(finalraw2[c("Braincode","Regioncode", "Window")],sep="_"))
finalraw2[(1:3)] <- NULL
finalraw2<-finalraw2 %>% select(which(colMedians(as.matrix(finalraw2))>2))
library(tibble)
finalraw2 <- rownames_to_column(finalraw2) 
finalraw2 <- finalraw2 %>% 
  separate(rowname, c("Braincode","Regioncode", "Window")) #Input to PC
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
colnames(df_raw2) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Ifan", "Infan/Child", "Child", "Adolescence", "Adult")
#PLOT
levels(colnames(df_raw2)) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Inf", "Inf/Child", "Child", "Adolescence", "Adult")
##write.csv(df_raw2, file="median_filtered_rawSestan.csv", row.names = FALSE)
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
##ggsave(file="filtered_rawSestan_log2_median.pdf", dfraw2_pl, width = 11.69, height = 8.27, units = "in")
```

#Akey in Sestan data
```{r}
mRNAseqData=read.table("~/tmp_psychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
modsb1= mRNAseqData %>% 
  separate(Geneid,c("EnsemblID","Genename"),extra="merge")
modsb1$EnsemblID<-NULL
#Filtering for Akey alone
akeySestan1 <- modsb1 %>% filter(modsb1$Genename %in% results$hgnc_symbol)
akeySestan1=t(akeySestan1)

#As dataframe
akeySestan1=as.data.frame(akeySestan1)

colnames(akeySestan1) <- as.matrix(unlist(akeySestan1[1,]))
akeySestan1 <- akeySestan1[-1, ]

akeySestan1 <- cbind(info = rownames(akeySestan1), akeySestan1)
rownames(akeySestan1) <- 1:nrow(akeySestan1)

#duplicated columns - issue in raw data. Here: 
colnames(akeySestan1)[duplicated(colnames(akeySestan1))] #1
akeySestan1 <- akeySestan1[, !duplicated(colnames(akeySestan1))]

akeySestan1=akeySestan1 %>% 
  separate(info, c("Braincode","Regioncode"))

#Normality
#normcheck <- akeypeySestan %>% pivot_longer(!("Braincode"|"Regioncode"), names_to = "Genes", values_to = "RPKM")
#normcheck$RPKM <- as.numeric(normcheck$RPKM)
#ggqqplot(normcheck, "RPKM", facet.by = "Regioncode")

#Transformation
logakeySestan1 <- akeySestan1

cols.num <- colnames(logakeySestan1[3:ncol(logakeySestan1)])

logakeySestan1[,cols.num] <- lapply(logakeySestan1[cols.num],as.character)
logakeySestan1[,cols.num] <- lapply(logakeySestan1[cols.num],as.numeric)

logakeySestan1[3:ncol(logakeySestan1)] <- log2(logakeySestan1[3:ncol(logakeySestan1)]+1)

#Normality after transformation
#lognormcheck <- logakeypeySestan %>% pivot_longer(!("Braincode"|"Regioncode"), names_to = "Genes", values_to = "RPKM")
#lognormcheck$RPKM <- as.numeric(lognormcheck$RPKM)
#ggqqplot(lognormcheck, "RPKM", facet.by = "Regioncode")


#Brining the metadata of the database
library(readxl)
metadatamRNAseq=read_xlsx("~/tmp_psychENCODE/mRNA-seq_QC.xlsx",skip = 3)


modMetadatamRNAseq=na.omit(metadatamRNAseq)
modMetadatamRNAseq=modMetadatamRNAseq %>% select(1:3)
modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)

#With log
finalakeySestan1=merge(modMetadatamRNAseq,logakeySestan1,by=c("Braincode", "Regioncode")) #input
 
# #write.csv(finalakeySestan1,"logakey.csv")

#Without log:
#finalakeySestan1=merge(modMetadatamRNAseq,akeypeySestan,by=c("Braincode", "Regioncode"))

finalakeySestan1$Braincode <- NULL
finalakeySestan1$Label<- paste(finalakeySestan1$Regioncode, finalakeySestan1$Window, sep="_")
finalakeySestan1$Regioncode <- NULL
finalakeySestan1$Window <- NULL

finalakeySestan1 <- as_tibble(finalakeySestan1)

finalakeySestan1_w_mean <- finalakeySestan1 %>%
    group_by(Label) %>%
    dplyr::summarise_all(median, na.rm=TRUE)

#For plotting
dfakeySestan1 <- data.frame(Structure=finalakeySestan1_w_mean$Label, Means=rowMeans(finalakeySestan1_w_mean[,-1]))

dfakeySestan1 <- dfakeySestan1 %>% 
  separate(Structure, c("Structure","Window"))
preli1 <- pivot_wider(dfakeySestan1, names_from = Window, values_from = Means)
preli1[,10] <- NULL
preli1 <- preli1[complete.cases(preli1), ]
colnames(preli1) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Ifan", "Infan/Child", "Child", "Adolescence", "Adult")
#PLOT
levels(colnames(preli1)) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Inf", "Inf/Child", "Child", "Adolescence", "Adult")

#write.csv(preli1, file="~/raul_tesina/1.data/median_tables/median_Akey_Sestan.csv", row.names = FALSE)

a<- ggparcoord(preli1,
              columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts- VFC (green) & AMY (black)")+scale_color_manual(values = c( "#ABABAB", "#000000", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#238b45"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
b <- ggparcoord(preli1,
                 columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
c<-ggparcoord(preli1,
               columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("expression")
a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)
# 
plak_sestan <-grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
#ggsave(file="~/raul_tesina/2.plots/Sestan_AkeyPey_log2_median/Sestan_GenesAkey_log2_median.pdf", plak_sestan, width = 11.69, height = 8.27, units = "in")

```

#AkeyPey in Sestan data
```{r}
mRNAseqData=read.table("~/tmp_psychENCODE/mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
modsb1= mRNAseqData %>% 
  separate(Geneid,c("EnsemblID","Genename"),extra="merge")
modsb1$EnsemblID<-NULL
#Filtering for Akey alone
akeypeySestan <- modsb1 %>% filter(modsb1$Genename %in% both$hgnc_symbol)
akeypeySestan=t(akeypeySestan)



#As dataframe
akeypeySestan=as.data.frame(akeypeySestan)

colnames(akeypeySestan) <- as.matrix(unlist(akeypeySestan[1,]))
akeypeySestan <- akeypeySestan[-1, ]

akeypeySestan <- cbind(info = rownames(akeypeySestan), akeypeySestan)
rownames(akeypeySestan) <- 1:nrow(akeypeySestan)

#duplicated columns - issue in raw data. Here: 
colnames(akeypeySestan)[duplicated(colnames(akeypeySestan))] #0
#akeypeySestan <- akeypeySestan[, !duplicated(colnames(akeypeySestan))]

akeypeySestan=akeypeySestan %>% 
  separate(info, c("Braincode","Regioncode"))

# #Normality
# normcheck <- akeypeySestan %>% pivot_longer(!("Braincode"|"Regioncode"), names_to = "Genes", values_to = "RPKM")
# normcheck$RPKM <- as.character(normcheck$RPKM)
# normcheck$RPKM <- as.numeric(normcheck$RPKM)
# 
# norm_1=merge(modMetadatamRNAseq,normcheck,by=c("Braincode", "Regioncode"))
# 
# norm_1 <- as_tibble(norm_1)
# 
# ggqqplot(norm_1, "RPKM", facet.by = "Window")
# 
# #Transformation
logakeypeysestan1 <- akeypeySestan
# 
cols.num <- colnames(logakeypeysestan1[3:ncol(logakeypeysestan1)])
# 
logakeypeysestan1[,cols.num] <- lapply(logakeypeysestan1[cols.num],as.character)
logakeypeysestan1[,cols.num] <- lapply(logakeypeysestan1[cols.num],as.numeric)
# 
logakeypeysestan1[3:ncol(logakeypeysestan1)] <- log2(logakeypeysestan1[3:ncol(logakeypeysestan1)]+1)
#ta <- logakeypeysestan1 %>% filter(logakeypeysestan1[3:ncol(logakeypeysestan1)]>1)

# 
# #Normality after transformation
# lognormcheck <- logakeypeysestan1 %>% pivot_longer(!("Braincode"|"Regioncode"), names_to = "Genes", values_to = "RPKM")
# lognormcheck$RPKM <- as.character(lognormcheck$RPKM)
# lognormcheck$RPKM <- as.numeric(lognormcheck$RPKM)
# lognorm_1=merge(modMetadatamRNAseq,lognormcheck,by=c("Braincode", "Regioncode"))
# # 
# lognorm_1 <- as_tibble(lognorm_1)
# # 
# ggqqplot(lognorm_1, "RPKM", facet.by = "Window")
# 

#Brining the metadata of the database
library(readxl)
metadatamRNAseq=read_xlsx("~/tmp_psychENCODE/mRNA-seq_QC.xlsx",skip = 3)


modMetadatamRNAseq=na.omit(metadatamRNAseq)
modMetadatamRNAseq=modMetadatamRNAseq %>% select(1:3)
modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)

#With log
finalakeypeySestan=merge(modMetadatamRNAseq,logakeypeysestan1,by=c("Braincode", "Regioncode"))
# #write.csv(finalakeypeySestan,"logakeypey.csv")
#Without log:
#finalakeypeySestan=merge(modMetadatamRNAseq,akeypeySestan,by=c("Braincode", "Regioncode"))

finalakeypeySestan$Braincode <- NULL
finalakeypeySestan$Label<- paste(finalakeypeySestan$Regioncode, finalakeypeySestan$Window, sep="_")
finalakeypeySestan$Regioncode <- NULL
finalakeypeySestan$Window <- NULL

finalakeypeySestan <- as_tibble(finalakeypeySestan)

finalakeypeySestan_w_mean <- finalakeypeySestan %>%
    group_by(Label) %>%
    dplyr::summarise_all(median, na.rm=TRUE)


#For plotting
dfakeypeySestan <- data.frame(Structure=finalakeypeySestan_w_mean$Label, Means=rowMeans(finalakeypeySestan_w_mean[,-1]))

dfakeypeySestan <- dfakeypeySestan %>% 
  separate(Structure, c("Structure","Window"))
preli2 <- pivot_wider(dfakeypeySestan, names_from = Window, values_from = Means)
preli2[,10] <- NULL
preli2 <- preli2[complete.cases(preli2), ]
colnames(preli2) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Ifan", "Infan/Child", "Child", "Adolescence", "Adult")
#PLOT
levels(colnames(preli2)) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Inf", "Inf/Child", "Child", "Adolescence", "Adult")

#write.csv(preli2, file="~/raul_tesina/1.data/median_tables/median_AkeyPey_Sestan.csv", row.names = FALSE)

a<- ggparcoord(preli2,
              columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts and Pos Sel- VFC (green) & AMY (black)")+scale_color_manual(values = c( "#ABABAB", "#000000", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#238b45"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
b <- ggparcoord(preli2,
                 columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Striatum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#FF0000", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10),legend.position = "none")+xlab("")+ylab("expression")
c<-ggparcoord(preli2,
               columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Cerebellum")+scale_color_manual(values = c( "#ABABAB", "#ABABAB", "#0000FF", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB",  "#ABABAB", "#ABABAB", "#ABABAB","#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB", "#ABABAB"))+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("expression")
a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a1, arrangeGrob(b1, c1), ncol = 2)
 
plakpey_sestan <-grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))

##ggsave(file="~/raul_tesina/2.plots/Sestan_AkeyPey_log2_median/Sestan_GenesAkeyPey_log2_median.pdf", plakpey_sestan, width = 11.69, height = 8.27, units = "in")
```

#Sestan - STATS
```{R}
#Friedman test per each stage independently.
##Post hoc via Conover
###AKEY & PEY
finalakeypeySestan_w_mean
fr <- finalakeypeySestan_w_mean %>% pivot_longer(!("Label"), names_to = "Genes", values_to = "expression")
fr <- fr %>% 
  separate(Label, c("Structure","Window"))
fr <- fr %>% filter(Window != 1)
fr <- fr %>% filter(Structure != "MSC")

ps <- lapply(split(fr, fr$Window), function(x) {friedman_test(expression ~ Structure | Genes, data=x)})

#pvalues <- numeric(0)
#for (i in (1:length(ps))) {pvalues[i] <- ps[[i]]$p}
#p.adjust(pvalues, method = "BH")

#Stages 3 & 4 & 7 
library(PMCMR)
fr3 <- fr %>% filter(fr$Window == 3)
a <-pivot_wider(fr3[,-2], names_from = Structure, values_from = expression)
a <- column_to_rownames(a, var = "Genes")
a <- as.matrix(a)
posthoc.friedman.conover.test(y=a, p.adjust="BH")

fr4 <- fr %>% filter(fr$Window == 4)
a <-pivot_wider(fr4[,-2], names_from = Structure, values_from = expression)
a <- column_to_rownames(a, var = "Genes")
a <- as.matrix(a)
posthoc.friedman.conover.test(y=a, p.adjust="BH")

fr7 <- fr %>% filter(fr$Window == 7)
a <-pivot_wider(fr7[,-2], names_from = Structure, values_from = expression)
a <- column_to_rownames(a, var = "Genes")
a <- as.matrix(a)
posthoc.friedman.conover.test(y=a, p.adjust="BH")

#Friedman test per each stage independently.
##Post hoc via Conover
###AKEY alone
finalakeySestan1_w_mean
fr0 <- finalakeySestan1_w_mean %>% pivot_longer(!("Label"), names_to = "Genes", values_to = "expression")
fr0 <- fr0 %>% 
  separate(Label, c("Structure","Window"))
fr0 <- fr0 %>% filter(Window != 1)
fr0 <- fr0 %>% filter(Structure != "MSC")

ps0 <- lapply(split(fr0, fr0$Window), function(x) {friedman_test(expression ~ Structure | Genes, data=x)})

#pvals0 <- numeric(0)
#for (i in (1:length(ps0))) {pvals0[i] <- ps0[[i]]$p}
#p.adjust(pvals0, method = "BH") #All stages report significant results
```

#Trajectory plots1 - ABA
```{r}
final_merge0 #all
final_merge #akey
final_mergelog #akeypey

final_merge0$dataset <- c("raw")
final_merge$dataset <- c("akey")
final_mergelog$dataset <- c("akeypey")

tot_pl <-rbind(final_merge, final_merge0, final_mergelog)
tot_pl <- as_tibble(tot_pl)
tot_pl <- tot_pl %>% mutate(dataset=as.character(dataset))

levels(tot_pl$dataset) <-  c("raw", "akeypey", "akey")

n <-ggparcoord(tot_pl,
columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Structures ABA")+theme(plot.title = element_text(size=10), legend.position = "none")+xlab("")+ylab("expression")



n <-ggparcoord(tot_pl,
columns = 2:6, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Structures ABA", mapping=aes(color=as.factor(dataset)))+theme(plot.title = element_text(size=10))+xlab("")+ylab("expression")+labs(color="Dataset")


n + facet_wrap(~Structure)

an <- n + facet_wrap(~Structure) 
##ggsave(file="~/raul_tesina/2.plots/ABAData_AkeyPeyRac_log2/ABA_temporal_Structures.pdf", an, width = 11.69, height = 8.27, units = "in")
```

#Trajectory plots2- Sestan
```{r}
preli1 <- read.csv(file="~/raul_tesina/1.data/median_tables/median_Akey_Sestan.csv")
preli1 <- as_tibble(preli1) #AKEY

preli2 <- read.csv(file="~/raul_tesina/1.data/median_tables/median_AkeyPey_Sestan.csv")
preli2 <- as_tibble(preli2) #AKEYpey

preli3 <- read.csv(file="~/raul_tesina/1.data/median_tables/median_filtered_rawSestan.csv")
preli3 <- as_tibble(preli3) #rawfiltered

colnames(preli1) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Ifan", "Infan/Child", "Child", "Adolescence", "Adult")
#PLOT
levels(colnames(preli1)) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Inf", "Inf/Child", "Child", "Adolescence", "Adult")

colnames(preli2) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Ifan", "Infan/Child", "Child", "Adolescence", "Adult")
#PLOT
levels(colnames(preli2)) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Inf", "Inf/Child", "Child", "Adolescence", "Adult")

colnames(preli3) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Ifan", "Infan/Child", "Child", "Adolescence", "Adult")
#PLOT
levels(colnames(preli3)) <- c("Structure", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Inf", "Inf/Child", "Child", "Adolescence", "Adult")

preli1$dataset <- c("Deserts")
preli2$dataset <- c("Deserts and Pos. Sel.")
preli3$dataset <- c("Global profile")

tot_pl <-rbind(preli1, preli2, preli3)
tot_pl <- as_tibble(tot_pl)
tot_pl <- tot_pl %>% mutate(dataset=as.character(dataset))

levels(tot_pl$dataset) <-  c("Deserts", "Deserts and Pos. Sel.", "Global profile")


n <-ggparcoord(tot_pl,
columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax", mapping=aes(color=as.factor(dataset)))+theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 45,  hjust = 1))+xlab("")+ylab("expression")+labs(color="Dataset")

n + facet_wrap(~Structure)+scale_color_discrete(name="Dataset",labels=unique(tot_pl$dataset))

an <- n + facet_wrap(~Structure)+scale_color_discrete(name="Dataset",labels=unique(tot_pl$dataset))
#ggsave(file="~/raul_tesina/2.plots/Sestan_raw_filtered_median2_trajectories/Sestan_temporal_Structures_filteredRawmedian2.pdf", an, width = 11.69, height = 8.27, units = "in")
```


#Preparing Figure X1
```{r}
traj_fig <- tot_pl %>% filter(Structure %in% "DFC" | Structure %in% "M1C" | Structure %in% "S1C" | Structure %in% "CBC") 

figtop1 <- ggparcoord(traj_fig,
           columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax", mapping=aes(color=as.factor(dataset)))+theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 45,  hjust = 1))+xlab("")+ylab("expression")+labs(color="Dataset")+ facet_wrap(~Structure, ncol = 4)+scale_color_discrete(name="Dataset",labels=unique(tot_pl$dataset))

#Akey
vpreli1 <- preli1
vpreli1$Structure <- str_replace(vpreli1$Structure, ("A1C|DFC|IPC|ITC|M1C|MFC|OFC|S1C| STC|V1C|VFC"), "NCX")

test <- vpreli1 %>%
    group_by(Structure) %>%
    dplyr::summarise_all(mean, na.rm=TRUE)
figtop2<-ggparcoord(test,
           columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="B")+theme(plot.title = element_text(size=10),legend.position = "none", axis.text.x = element_text(angle = 45,  hjust = 1))+xlab("")+ylab("expression")+geom_line(size=1.7, alpha=0.7)


# figtop2 <- ggparcoord(vpreli1,
#               columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Deserts of introgression")+theme(plot.title = element_text(size=10),legend.position = "none", axis.text.x = element_text(angle = 45,  hjust = 1))+xlab("")+ylab("expression")
#AKeyPey
vpreli2 <- preli2
vpreli2$Structure <- str_replace(vpreli2$Structure, ("A1C|DFC|IPC|ITC|M1C|MFC|OFC|S1C| STC|V1C|VFC"), "NCX")

test2 <- vpreli2 %>%
    group_by(Structure) %>%
    dplyr::summarise_all(mean, na.rm=TRUE)
figtop3<-ggparcoord(test2,
           columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="C")+theme(plot.title = element_text(size=10),legend.position = "right", axis.text.x = element_text(angle = 45,  hjust = 1))+xlab("")+ylab("expression")+geom_line(size=1.7, alpha=0.7)

# figtop3 <- ggparcoord(preli2,
#               columns = 2:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Deserts and Positively selected regions")+theme(plot.title = element_text(size=10),legend.position = "right", axis.text.x = element_text(angle = 45,  hjust = 1))+xlab("")+ylab("expression")

lay <- rbind(c(1,1,1,1,1,1),
             c(2,2,2,3,3,3), 
             c(2,2,2,3,3,3))

grid.arrange(figtop1, figtop2, figtop3, layout_matrix = lay)
```

#Preparing Figure X2
```{r}
wilcoxAkeyDist <- read.csv("~/tmp_wilcoxtests/wilcoxAkeyDist.csv")




```



#Any set of genes
```{r}
ap <- finalakeySestan1 %>% select(both$hgnc_symbol|Label)
tes2 <- melt(ap, id=c("Label"))
tes2_mean <- tes2 %>%
  group_by(Label,variable) %>%
  dplyr::summarise_all(median, na.rm=TRUE)

tes2_mean <- tes2_mean %>% 
  separate(Label, c("Structure","Window"))
prelites2 <- pivot_wider(tes2_mean, names_from = Window, values_from = value)

#
prelites2[,11] <- NULL
prelites2 <- prelites2[complete.cases(prelites2), ]

#levels(prelites2$Genename) <- unique(prelites2$Genename)

colnames(prelites2) <- c("Structure", "Genename", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Ifan", "Infan/Child", "Child", "Adolescence", "Adult")
#PLOT
levels(colnames(prelites2)) <- c("Structure", "Genename", "Fetal_1", "Fetal_2", "Fetal_3", "Birth/Inf", "Inf/Child", "Child", "Adolescence", "Adult")

n <-ggparcoord(prelites2,
columns = 3:9, groupColumn = 2, showPoints = TRUE, scale = "globalminmax", mapping=aes(color=factor(Structure)))+xlab("")+ylab("expression")+theme(axis.text.x = element_text(angle = 45,hjust = 1))
j1<-n + facet_wrap(~Genename)+scale_color_discrete(name="Structure",labels=unique(prelites2$Structure))
#ggsave(file="~/raul_tesina/2.plots/Sestan_AkeyPey_log2_median//GenesAkeyPey_log2_median_perGenes_1.pdf", j1, width = 11.69, height = 8.27, units = "in")

n1 <-ggparcoord(prelites2,
columns = 3:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts and Pos Sel", mapping=aes(color=factor(Genename)))+xlab("")+ylab("expression")
j2 <- n1 + facet_wrap(~Structure)+scale_color_discrete(name="Genes",labels=unique(prelites2$Genename))+theme(axis.text.x = element_text(angle = 45,hjust = 1))
#ggsave(file="~/raul_tesina/2.plots/Sestan_AkeyPey_log2_median/GenesAkeyPey_log2_median_perGenes_2.pdf", j2, width = 11.69, height = 8.27, units = "in")

#Say CBC e.g.
# t1 <- filter(prelites2, Structure %in% c("CBC")) %>% ggparcoord(prelites2,
# columns = 3:9, groupColumn = 1, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts and Pos Sel", mapping=aes(color=factor(Genename)))+xlab("")+ylab("expression")
# 
# t1 + facet_wrap(~Structure)+scale_color_discrete(name="Genes",labels=unique(prelites2$Genename))
# 
# t2 <- filter(prelites2, Structure %in% c("CBC")) %>% ggparcoord(prelites2,
# columns = 3:9, groupColumn = 2, showPoints = TRUE, scale = "globalminmax",title="Genes in Deserts and Pos Sel", mapping=aes(color=factor(Structure)))+xlab("")+ylab("expression")
# t2 + facet_wrap(~Genename)+scale_color_discrete(name="Structure",labels=c("CBC"))
```
