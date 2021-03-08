#ABA PeyNotAkey
#Preparing ABAdata for PCA and distances analysis

##Loading dataset
data("dataset_5_stages")
resultsAkey<-read.csv("results_akey.csv")
resultsAkeyPey<-read.csv("both_pey_akey_genes_pos.csv")
#Selecting genes ID and structures present in dataset
id <- unique(dataset_5_stages$ensembl_gene_id)
st <- unique(dataset_5_stages$structure)
st_allen <- paste("Allen",st, sep=":") 
#Expression data for all structures and genes
ab <- get_expression(structure_ids=st_allen, gene_ids = id, dataset='5_stages')
abStatic<-copy(ab)
#Converting data to a dataframe with useful format for later
list1 = vector(mode="list")
for (r in 1:length(ab)){
  ab[[r]] <- t(ab[[r]]) #transpose
  list1 <- get_name(colnames(ab[[r]])) #change Allen:XXXX to e.g. M1C_primary motor cortex, etc
  colnames(ab[[r]]) <- list1
  ab[[r]] <- as.data.frame(ab[[r]])
}

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


peyNotakeyABA = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    peyNotakeyABA[[i]][[h]] <- ab1[[i]][[h]][ab1[[i]][[h]][[2]] %in% resultsPeynotAkey$hgnc_symbol,] #in Akey
    peyNotakeyABA[[i]][[h]] <-  peyNotakeyABA[[i]][[h]][!(is.na(peyNotakeyABA[[i]][[h]][[2]]) | peyNotakeyABA[[i]][[h]][[2]]==""), ]
    # akeyABA[[i]][[h]][1] <- log2(akeyABA[[i]][[h]][1]+1)
    #Cleaning
  }
}

peyNotakeyABAwind<-list()
peyNotakeyABAwind[["prenatal"]]<-as.data.frame(peyNotakeyABA[[1]])
peyNotakeyABAwind[["infant"]]<-as.data.frame(peyNotakeyABA[[2]])
peyNotakeyABAwind[["child"]]<-as.data.frame(peyNotakeyABA[[3]])
peyNotakeyABAwind[["adolsecent"]]<-as.data.frame(peyNotakeyABA[[4]])
peyNotakeyABAwind[["adult"]]<-as.data.frame(peyNotakeyABA[[5]])
#For each window handling
for(i in 1:5){
  peyNotakeyABAdf<-peyNotakeyABAwind[[i]]
  peyNotakeyABAdfTemp <- peyNotakeyABAdf %>% select(-contains("gene_name"))
  peyNotakeyABAdfTemp<-as.data.frame(t(peyNotakeyABAdfTemp))
  peyNotakeyABAdfTemp <- cbind(Regioncode = rownames(peyNotakeyABAdfTemp), peyNotakeyABAdfTemp)
  rownames(peyNotakeyABAdfTemp) <- 1:nrow(peyNotakeyABAdfTemp)
  peyNotakeyABAdfTemp<-peyNotakeyABAdfTemp%>%separate(Regioncode,c("Regioncode","Extra"),sep = "_")
  peyNotakeyABAdfTemp$Extra<-NULL
  peyNotakeyABAwind[[i]]<-peyNotakeyABAdfTemp
}



#PeynotAkey in Sestan data
mRNAseqData=read.table("mRNA-seq_hg38.gencode21.wholeGene.geneComposite.STAR.nochrM.gene.RPKM.normalized.CQNCombat.txt",sep="\t",header=TRUE)
modsb1= mRNAseqData %>% 
  separate(Geneid,c("EnsemblID","Genename"),extra="merge")
modsb1$EnsemblID<-NULL
#Filtering for Akey alone
# PeynotAkeySestan <- modsb1 %>% filter(modsb1$Genename %in% results$hgnc_symbol)
PeynotAkeySestan <- modsb1 %>% filter(modsb1$Genename %in% resultsPeynotAkey$hgnc_symbol)
resultsPeynotAkey
PeynotAkeySestan=t(PeynotAkeySestan)

PeynotAkeySestan=na.omit(PeynotAkeySestan)

#As dataframe
PeynotAkeySestan=as.data.frame(PeynotAkeySestan)

colnames(PeynotAkeySestan) <- as.matrix(unlist(PeynotAkeySestan[1,]))
PeynotAkeySestan <- PeynotAkeySestan[-1, ]

PeynotAkeySestan <- cbind(info = rownames(PeynotAkeySestan), PeynotAkeySestan)
rownames(PeynotAkeySestan) <- 1:nrow(PeynotAkeySestan)

#duplicated columns - issue in raw data. Here: 
colnames(PeynotAkeySestan)[duplicated(colnames(PeynotAkeySestan))] #0
PeynotAkeySestan <- PeynotAkeySestan[, !duplicated(colnames(PeynotAkeySestan))]

PeynotAkeySestan=PeynotAkeySestan %>% 
  separate(info, c("Braincode","Regioncode"))

# #Normality
# normcheck <- PeynotAkeySestan %>% pivot_longer(!("Braincode"|"Regioncode"), names_to = "Genes", values_to = "RPKM")
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
logPeynotAkeysestan1 <- PeynotAkeySestan
# 
cols.num <- colnames(logPeynotAkeysestan1[3:ncol(logPeynotAkeysestan1)])
# 
logPeynotAkeysestan1[,cols.num] <- lapply(logPeynotAkeysestan1[cols.num],as.character)
logPeynotAkeysestan1[,cols.num] <- lapply(logPeynotAkeysestan1[cols.num],as.numeric)
# 
logPeynotAkeysestan1[3:ncol(logPeynotAkeysestan1)] <- log2(logPeynotAkeysestan1[3:ncol(logPeynotAkeysestan1)]+1)
#ta <- logPeynotAkeysestan1 %>% filter(logPeynotAkeysestan1[3:ncol(logPeynotAkeysestan1)]>1)
write.csv(logPeynotAkeysestan1,"peyNotAkeySeastanLogNorm.csv")

#PCA willcox and boxplots
```{r}
akeyPeyNorm=read.csv("logakeypey.csv")
akeyNorm=read.csv("logakey.csv")
PeynotAkeySestan=read.csv("peyNotAkeySeastanLogNorm.csv")

logakeypeyPCA<-akeyPeyNorm
logakeypeyPCA$X<-NULL
logakeypeyPCA$Braincode<-NULL

logakeyPCA<-akeyNorm
logakeyPCA$X<-NULL
logakeyPCA$Braincode<-NULL

logPeynotAkeySestan<-PeynotAkeySestan
logPeynotAkeySestan$X<-NULL
logPeynotAkeySestan$Braincode<-NULL

library(readxl)
metadatamRNAseq=read_xlsx("mRNA-seq_QC.xlsx",skip = 3)


modMetadatamRNAseq=na.omit(metadatamRNAseq)
modMetadatamRNAseq=modMetadatamRNAseq %>% select(1:3)
modMetadatamRNAseq=as.data.frame(modMetadatamRNAseq)
lognorm_1=merge(modMetadatamRNAseq,PeynotAkeySestan,by=c("Braincode", "Regioncode"))
# 
lognorm_1 <- as_tibble(lognorm_1)
lognorm_1<-na.omit(lognorm_1)
# 
ggqqplot(lognorm_1, "RPKM", facet.by = "Window")

# write.csv(lognorm_1,"peyNotAkeySeastanLogNorm.csv")
logPeynotAkeySestan<-na.omit(lognorm_1)
lognorm_1[,4:ncol(lognorm_1)]<-lapply(lognorm_1[,4:ncol(lognorm_1)],as.numeric)
logPeynotAkeySestan<-lognorm_1
logPeynotAkeySestan$X<-NULL
logPeynotAkeySestan$Braincode<-NULL
structDistWind<-list()
for (i in 1:5){
  
  #CHANGE:logakeypeyPCA / logakeyPCA / logPeynotAkeySestan
  # windowPCA<-logPeynotAkeySestan %>% filter(Window==i)
  
  #CHANGE: akeyABAwind / akeyPeyABAwind / peyNotAkeyABAwind / wholeABAwind
  windowPCA<-peyNotakeyABAwind[[i]]
  
  windowPCA$Window<-NULL
  
  
  brainTopStruct<-read.csv("brainRegionCorresp.csv",sep=":")
  windowPCASt=merge(brainTopStruct,windowPCA,by="Regioncode")
  #Cleaning columns equal to 0
  
  windowPCAStTemp<-windowPCASt[,-1][,-1][,-1][,-1][,-(which(colSums(windowPCASt[,-1][,-1][,-1][,-1])==0))]
  if(length(windowPCAStTemp)==0){
    windowPCAStTemp<-windowPCASt[,-1][,-1][,-1][,-1]
  }
  pca_res <- prcomp(windowPCAStTemp, scale. = TRUE)
  # PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$topStructure) 
  PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$brainAreas) #with neocortex
  structDist<-hash()
  for (structure in unique(PCi$topStructure)){
    
    dfStructure<-PCi %>% filter(topStructure==structure)
    dfOtherStructures<-PCi %>% filter(topStructure!=structure)
    distancesStruct<-c()
    sdStruct<-c()
    for (row in 1:nrow(dfStructure)) {
      xstruct<-dfStructure[row,]$PC1
      ystruct<-dfOtherStructures[row,]$PC2
      nrowsOthers<-nrow(dfOtherStructures)
      listDist<-c()
      for (rowOth in 1:nrowsOthers){
        xOtherstruct<-dfOtherStructures[rowOth,]$PC1
        yOtherstruct<-dfOtherStructures[rowOth,]$PC2
        distance<-dist(matrix(c(xstruct,ystruct,xOtherstruct,yOtherstruct),nrow=2,ncol=2),diag=TRUE)
        #print(c("Distance:",distance))
        listDist<-append(listDist,distance[1])
        nrowsOthers
      }
      distancesStruct<-append(distancesStruct,listDist)
      sdStruct<-append(sdStruct,sd(listDist))
      print(paste("Distance ",structure,": list->",head(listDist),"; sd->",sd(listDist)))
    }
    structDist[[structure]]<-c(distancesStruct,mean(sdStruct))
  }
  
  
  structDistWind[[i]]<-structDist
}

# Final var structDistWind
rm(pairwiseWilcox)

#Initialize dataframe
#ABA ncol=30 and 1:5 //SESTAN ncol=30 and 1:5
wilcoxTests<-data.frame(matrix(ncol=30,nrow = 6))
wilcoxTestsCol<-c()
for (wind in 1:5){
  for (str in keys(structDist)){
    wilcoxTestsCol<-append(wilcoxTestsCol,paste(str,"_",wind))
  }
}
colnames(wilcoxTests)<-wilcoxTestsCol
rownames(wilcoxTests)<-keys(structDist)
pairwiseWilcox<-list()
for (i in 1:5){
  print(i)
  structDist<-structDistWind[[i]]
  structs<-c()
  valuesDist<-c()
  for(str in keys(structDist)){
    structs<-append(structs,rep(str,length(values(structDist[str]))))
    valuesDist<-append(valuesDist,values(structDist[str]))
    #Other way to calculate manually multiple sequentially
    # for(str2 in keys(structDist)){
    #   if(str==str2){
    #     wilcoxTests[str2,paste(str,"_",i)]<-NaN
    #   }else{
    #     wilcoxTests[str2,paste(str,"_",i)]<-wilcox.test(values(structDist[str]),values(structDist[str2]))$p.value
    #   }
    # }
    
  }
  pairwiseWilcox[[correspStage[[i]]]]<-as.data.frame(pairwise.wilcox.test(valuesDist,structs, p.adj = "bonf")$p.value)
}
valuesDist
structs
library(xlsx)
install.packages("xlsx")
#wilcoxAkeyDist / wilcoxAkeyPeyDist
#writexl::write_xlsx(pairwiseWilcox,"wilcoxPeyNotAkeyDist.xlsx")
plot(pairwiseWilcox[[2]])
#colMeans(wilcoxTests, na.rm = TRUE)

#Visualizing p-values wilcox

#First handling the data
wilcoxTests<-data.frame(matrix(ncol=30,nrow = 6))
wilcoxTestsCol<-c()
for (wind in 1:5){
  for (str in keys(structDist)){
    wilcoxTestsCol<-append(wilcoxTestsCol,paste(str,"_",wind))
  }
}
colnames(wilcoxTests)<-wilcoxTestsCol
rownames(wilcoxTests)<-keys(structDist)
for (i in 1:5){
  #Restructuring the data
  #Introducing to the dataframe
  
  isSaved<-list()
  pairwiseData<-pairwiseWilcox[[correspStage[[i]]]]
  for(str in keys(structDist)){
    for(str2 in keys(structDist)){
      print(paste(str,"-->",str2,": ",pairwiseData[str,str2],is.null(pairwiseData[str,str2]),is.na(pairwiseData[str,str2])))
      #isSaved[[paste(str,str2)]-->isSaved[[paste(str2,str)] for just half of the tables
      if((is.null(isSaved[[paste(str,str2)]]))&&(str!=str2)){
        if(!is.null(pairwiseData[str,str2])){
          if(!is.na(pairwiseData[str,str2])){
            wilcoxTests[str2,paste(str,"_",i)]<-pairwiseData[str,str2]
          }else{
            if(!is.null(pairwiseData[str2,str])){
              if(!is.na(pairwiseData[str2,str])){
                wilcoxTests[str2,paste(str,"_",i)]<-pairwiseData[str2,str]
              }else{
                wilcoxTests[str2,paste(str,"_",i)]<-NA
              }
            }else{
              wilcoxTests[str2,paste(str,"_",i)]<-NA
            }
          }
        } else if(!is.null(pairwiseData[str2,str])){
          if(!is.na(pairwiseData[str2,str])){
            wilcoxTests[str2,paste(str,"_",i)]<-pairwiseData[str2,str]
          }else{
            wilcoxTests[str2,paste(str,"_",i)]<-NA
          }
        }else{
          wilcoxTests[str2,paste(str,"_",i)]<-NA
        }
        
      }else{
        wilcoxTests[str2,paste(str,"_",i)]<-NA
      }
      isSaved[[paste(str,str2)]]<-TRUE
      
    }
  }
}

willcoxPvaluesAVG<-as.data.frame(t(colMeans(wilcoxTests,na.rm=TRUE)))

# wilcoxAkeyDist.csv/wilcoxAkeyPeyDist.csv
write.csv(wilcoxTests,"wilcoxABAPeynotAkeyDist.csv")

willCoxPvals<-as.data.frame(list("brainRegion","window","pvalAVG"))
colnames(willCoxPvals)<-c("brainRegion","window","pvalAVG")
willCoxPvals$pvalAVG<-as.numeric(willCoxPvals$pvalAVG)
willCoxPvals<-willCoxPvals[-1,]
for (col in colnames(willcoxPvaluesAVG)){
  #change the window number
  willCoxPvals<-willCoxPvals %>% add_row(brainRegion=strsplit(col," _ ")[[1]][1],window=strsplit(col," _ ")[[1]][2],pvalAVG=willcoxPvaluesAVG[col][[1]])
}
willCoxPvals$log2pvalAVG<-log2(willCoxPvals$pvalAVG)
plotWillCox<-ggplot(willCoxPvals, aes(x=window, y=log2pvalAVG, group=brainRegion)) +
  geom_line(aes(color=brainRegion))+
  geom_point(aes(color=brainRegion))
plotWillCox
ggsave(file="ABA_WillcoxPvalPlotDistances_PeynotAkey.pdf", width = 11.69, height = 8.27)

#Boxplots
library(ggplot2)
library(reshape2)

correspStage<-list()
# correspStage[[2]]<-"fetal1"
# correspStage[[3]]<-"fetal2"
# correspStage[[4]]<-"fetal3"
# correspStage[[5]]<-"Birth/Infan"
# correspStage[[6]]<-"Infan/Chil"
# correspStage[[7]]<-"Child"
# correspStage[[8]]<-"Adolescent"
# correspStage[[9]]<-"Adult"
correspStage[[1]]<-"Prenatal"
correspStage[[2]]<-"Infant"
correspStage[[3]]<-"Child"
correspStage[[4]]<-"Adolescent"
correspStage[[5]]<-"Adult"
library(reshape2)
boxplotsDist<-list()
for (i in 1:5){
  structDist<-structDistWind[[i]]
  valoresStruct<-list()
  for (structure in keys(structDist)){
    valoresStruct[[structure]]<-values(structDist[structure])
  }
  #In a melting pot
  valoresStructdf <- melt(valoresStruct)
  valoresStructdf$Var1<-NULL
  valoresStructdf$L1<-NULL
  
  # prepare a special xlab with the number of obs for each group
  my_xlab <- paste(levels(valoresStructdf$Var1),"\n(N=",table(valoresStructdf$Var1),")",sep="")
  colnames(valoresStructdf) <- c("names", "expression")
  # plot
  boxplotsDist[[i]]<-ggplot(valoresStructdf, aes(x=names, y=expression, fill=names)) +
    geom_boxplot(varwidth = TRUE) +
    theme(legend.position="none",axis.text.x = element_blank()) + xlab(correspStage[[i]])
  
}
#ABA
ggarrange(boxplotsDist[[1]],boxplotsDist[[2]],boxplotsDist[[3]],
          boxplotsDist[[4]],boxplotsDist[[5]],plotWillCox,
          common.legend = TRUE, legend = "right")
# Sestan
# ggarrange(boxplotsDist[[2]],boxplotsDist[[3]],
#           boxplotsDist[[4]],boxplotsDist[[5]],boxplotsDist[[6]],boxplotsDist[[7]],boxplotsDist[[8]],boxplotsDist[[9]],plotWillCox,
#           common.legend = TRUE, legend = "right")
ggsave(file="ABA_peyNotAkey_supplfig.pdf", width = 11.69, height = 8.27)


