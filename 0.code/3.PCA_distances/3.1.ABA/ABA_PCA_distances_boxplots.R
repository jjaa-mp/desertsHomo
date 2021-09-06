#Preparing ABAdata for PCA and distances analysis
library(tsne)
library(biomaRt)
library(Seurat)
library(ggplot2)
library(fviz_pca_ind)
library("factoextra")
library(ABAData)
library(ABAEnrichment)
library(data.table)
library(ggplot2)
library(dplyr)
library(GGally)
library(viridis)
library(factoextra)
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
library(cluster)
library(purrr)

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
##Pey coordinates (Ensembl 1-based):
pey_coords <- read.delim("desertsHomo/1.data/input_data/2020_pey_coords.bed", header=FALSE) #File with start 0-based
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


PeyNotAkey<- results_pey[!results_pey$hgnc_symbol %in% results$hgnc_symbol,]
PeyNotAkey <-  PeyNotAkey[!(is.na(PeyNotAkey$hgnc_symbol) | PeyNotAkey$hgnc_symbol==""), ] #Cleaning
##Loading dataset
data("dataset_5_stages")
resultsAkey<-results # read.csv("results_akey.csv")
resultsAkeyPey<-both # read.csv("both_pey_akey_genes_pos.csv")
resultsPeynotAkey<-PeyNotAkey # read.csv("both_pey_akey_genes_pos.csv")
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



akeyABA = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    akeyABA[[i]][[h]] <- ab1[[i]][[h]][ab1[[i]][[h]][[2]] %in% results$hgnc_symbol,] #in Akey
    akeyABA[[i]][[h]] <-  akeyABA[[i]][[h]][!(is.na(akeyABA[[i]][[h]][[2]]) | akeyABA[[i]][[h]][[2]]==""), ]
    # akeyABA[[i]][[h]][1] <- log2(akeyABA[[i]][[h]][1]+1)
    #Cleaning
  }
}

akeyABAwind<-list()
akeyABAwind[["prenatal"]]<-as.data.frame(akeyABA[[1]])
akeyABAwind[["infant"]]<-as.data.frame(akeyABA[[2]])
akeyABAwind[["child"]]<-as.data.frame(akeyABA[[3]])
akeyABAwind[["adolsecent"]]<-as.data.frame(akeyABA[[4]])
akeyABAwind[["adult"]]<-as.data.frame(akeyABA[[5]])
#For each window handling
for(i in 1:5){
  akeyABAdf<-akeyABAwind[[i]]
  akeyABAdfTemp <- akeyABAdf %>% select(-contains("gene_name"))
  akeyABAdfTemp<-as.data.frame(t(akeyABAdfTemp))
  akeyABAdfTemp <- cbind(Regioncode = rownames(akeyABAdfTemp), akeyABAdfTemp)
  rownames(akeyABAdfTemp) <- 1:nrow(akeyABAdfTemp)
  akeyABAdfTemp<-akeyABAdfTemp%>%separate(Regioncode,c("Regioncode","Extra"),sep = "_")
  akeyABAdfTemp$Extra<-NULL
  akeyABAwind[[i]]<-akeyABAdfTemp
}


wholeABA = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    wholeABA[[i]][[h]] <- ab1[[i]][[h]][ab1[[i]][[h]][[2]] %in% ab1[[i]][[h]][[2]],] #in Akey
    wholeABA[[i]][[h]] <-  wholeABA[[i]][[h]][!(is.na(wholeABA[[i]][[h]][[2]]) | wholeABA[[i]][[h]][[2]]==""), ]
    # akeyPeyABA[[i]][[h]][1] <- log2(akeyPeyABA[[i]][[h]][1]+1)
    #Cleaning
  }
}
wholeABAwind<-list()
wholeABAwind[["prenatal"]]<-as.data.frame(wholeABA[[1]])
wholeABAwind[["infant"]]<-as.data.frame(wholeABA[[2]])
wholeABAwind[["child"]]<-as.data.frame(wholeABA[[3]])
wholeABAwind[["adolsecent"]]<-as.data.frame(wholeABA[[4]])
wholeABAwind[["adult"]]<-as.data.frame(wholeABA[[5]])
#For each window handling
for(i in 1:5){
  wholeABAdf<-wholeABAwind[[i]]
  wholeABAdfTemp <- wholeABAdf %>% select(-contains("gene_name"))
  wholeABAdfTemp<-as.data.frame(t(wholeABAdfTemp))
  wholeABAdfTemp <- cbind(Regioncode = rownames(wholeABAdfTemp), wholeABAdfTemp)
  rownames(wholeABAdfTemp) <- 1:nrow(wholeABAdfTemp)
  wholeABAdfTemp<-wholeABAdfTemp%>%separate(Regioncode,c("Regioncode","Extra"),sep = "_")
  wholeABAdfTemp$Extra<-NULL
  wholeABAwind[[i]]<-wholeABAdfTemp
}
akeyPeyABA = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    akeyPeyABA[[i]][[h]] <- ab1[[i]][[h]][ab1[[i]][[h]][[2]] %in% both$hgnc_symbol,] #in Akey
    akeyPeyABA[[i]][[h]] <-  akeyPeyABA[[i]][[h]][!(is.na(akeyPeyABA[[i]][[h]][[2]]) | akeyPeyABA[[i]][[h]][[2]]==""), ]
    # akeyPeyABA[[i]][[h]][1] <- log2(akeyPeyABA[[i]][[h]][1]+1)
    #Cleaning
  }
}

akeyPeyABAwind<-list()
akeyPeyABAwind[["prenatal"]]<-as.data.frame(akeyPeyABA[[1]])
akeyPeyABAwind[["infant"]]<-as.data.frame(akeyPeyABA[[2]])
akeyPeyABAwind[["child"]]<-as.data.frame(akeyPeyABA[[3]])
akeyPeyABAwind[["adolsecent"]]<-as.data.frame(akeyPeyABA[[4]])
akeyPeyABAwind[["adult"]]<-as.data.frame(akeyPeyABA[[5]])
#For each window handling
for(i in 1:5){
  akeyPeyABAdf<-akeyPeyABAwind[[i]]
  akeyPeyABAdfTemp <- akeyPeyABAdf %>% select(-contains("gene_name"))
  akeyPeyABAdfTemp<-as.data.frame(t(akeyPeyABAdfTemp))
  akeyPeyABAdfTemp <- cbind(Regioncode = rownames(akeyPeyABAdfTemp), akeyPeyABAdfTemp)
  rownames(akeyPeyABAdfTemp) <- 1:nrow(akeyPeyABAdfTemp)
  akeyPeyABAdfTemp<-akeyPeyABAdfTemp%>%separate(Regioncode,c("Regioncode","Extra"),sep = "_")
  akeyPeyABAdfTemp$Extra<-NULL
  akeyPeyABAwind[[i]]<-akeyPeyABAdfTemp
}

peyNotakeyABA = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
for (i in 1:length(ab)){ #(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
  for (h in 1:length(names(ab[[i]]))){
    peyNotakeyABA[[i]][[h]] <- ab1[[i]][[h]][ab1[[i]][[h]][[2]] %in% PeyNotAkey$hgnc_symbol,] #in Akey
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

#PCA AkeyPey

#ABA
#akeyPeyABAwind here you can calculate the PCA for each window for the selected sets akeyABAwind/akeyPeyABAwind/peyNotAkeyABAwind/wholeABAwind


correspStage<-list()
correspStage[[1]]<-"Prenatal"
correspStage[[2]]<-"Infant"
correspStage[[3]]<-"Child"
correspStage[[4]]<-"Adolescent"
correspStage[[5]]<-"Adult"

jackSplot<-list()
for( i in 1:5){
  #CHANGE:akeyABAwind/akeyPeyABAwind/peyNotakeyABAwind/wholeABAwind
  windowPCA<-wholeABAwind[[i]] #example
  windowPCA$Window<-NULL
  brainTopStruct<-read.csv("desertsHomo/0.code/brainRegionCorresp.csv",sep=":")
  windowPCASt=merge(brainTopStruct,windowPCA,by="Regioncode")
  
  
  counts<-data.frame(t(windowPCASt[,-1][,-1][,-1][,-1]))
  meta<-data.frame(t(data.frame(t(windowPCASt[1:4]))))
  seurat_obj<-CreateSeuratObject(counts, project = "SeuratProject", assay = "RNA", meta.data = meta)
  #We have 16 features so we define the maximum number of PCS as 15
  numpcs<-15
  
  #Then at some point we have 12 items in one case, we define in it in 11 pcs
  if(nrow(counts)==12){
    numpcs<-11
  }
  all.genes <- rownames(counts)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  seurat_obj<-  FindVariableFeatures(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj,npcs =numpcs)
  seurat_obj <- JackStraw(seurat_obj, num.replicate = 100)
  seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:numpcs)
  jackSplot[[i]]<-JackStrawPlot(seurat_obj, dims = 1:numpcs)+ggtitle(correspStage[[i]])

}
correspStage[[i]]
tsne_res <- tsne(windowPCASt[,-1][,-1][,-1][,-1])
ev <- pca_res$sdev^2
pca_res.ve <- ev/sum(ev)
varPCA<-{par(mfrow = c(1,2), mar = c(4,5,3,1))
  plot(pca_res.ve,
       xlab = "Principal Component",
       ylab = "Proportion of Variance Explained",  
       type = 'b',
       main = paste('Scree plot,','Krule: ',toString(length(ev[ev>=1]))))
  
  plot(cumsum(pca_res.ve), 
       xlab = "Principal Component", 
       ylab = "Cumulative Proportion of\nVariance Explained", 
       type = 'b',
       main = paste('Scree plot,','Krule: ',toString(length(ev[ev>=1]))))
}
numPCA<-length(ev[ev>=1]) # Number of selected PCA by Kaiser’s rule
#With alternative sets

# PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$topStructure)
PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$brainAreas) #Neocortex

autoplot(pca_res, data = PCi, colour = 'Species', shape = FALSE,label.size=6)
tsnei<-data.frame(tsne_res,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$brainAreas) #Neocortex
tsnei
set.seed(123)


#Generating the plot of the pca for each window: w1,w2,w3,w4,w5
w5<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()#+
  #   theme(legend.position = "none")

ind.p <- fviz_pca_ind(pca_res, geom = "point", col.ind = PCi$topStructure)
plot2<-ggpubr::ggpar(ind.p,
              title = ,
              caption = "Source: factoextra",
              xlab = "PC1", ylab = "PC2",
              legend.title = "Species", legend.position = "top",
              ggtheme = theme_gray(), palette = "jco"
)
ggarrange(plot1,plot2)
#Joining plots
# ggarrange(w1, w2,w3,w4,w5, labels = c("prenatal", "infant","child","adolsecent","adult"),
#          common.legend = TRUE, legend = "right")

# ggsave(file="ABA_PCA_WindowsAkeyPeyNCXnoLog.pdf", width = 11.69, height = 8.27) #axample

#Euclidean distance between brain regions


correspStage<-list()
correspStage[[1]]<-"Prenatal"
correspStage[[2]]<-"Infant"
correspStage[[3]]<-"Child"
correspStage[[4]]<-"Adolescent"
correspStage[[5]]<-"Adult"

library(hash)
plotsKmeans<-list()
plotsPCA<-list()
varPCA<-list()
structDistWind<-list()
for (i in 1:5){
  
#CHANGE:akeyABAwind/akeyPeyABAwind/peyNotakeyABAwind/wholeABAwind
  
  windowPCA<-akeyABAwind[[i]]
  
  windowPCA$Window<-NULL
  
  
  brainTopStruct<-read.csv("desertsHomo/0.code/brainRegionCorresp.csv",sep=":")
  windowPCASt=merge(brainTopStruct,windowPCA,by="Regioncode")
  #Cleaning columns equal to 0
  windowPCAStTemp<-windowPCASt[,-1][,-1][,-1][,-1][,-(which(colSums(windowPCASt[,-1][,-1][,-1][,-1])==0))]
  if(length(windowPCAStTemp)==0){
    windowPCAStTemp<-windowPCASt[,-1][,-1][,-1][,-1]
  }
  
  pca_res <- prcomp(windowPCAStTemp, scale. = TRUE)
  #tsne_res <- tsne(windowPCAStTemp) # Alternative to PCA
  PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$brainAreas) 
  #PCi<-data.frame(tsne_res,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$brainAreas) #with neocortex
  
  #PCA plots
  ind.p <- fviz_pca_ind(pca_res, geom = "point", col.ind = PCi$topStructure)
  plotsPCA[[i]]<-ggpubr::ggpar(ind.p,
                       title = correspStage[[i]],
                       caption = "Source: factoextra",
                       xlab = "PC1", ylab = "PC2",
                       legend.title = "Species", legend.position = "top",
                       ggtheme = theme_gray(), palette = "jco"
  )
  #tsnei<-data.frame(tsne_res,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$brainAreas) #Neocortex

  ev <- pca_res$sdev^2
 
  # res.km <- kmeans(scale(windowPCAStTemp), 3, nstart = 25)
  # 
  # ind.coord <- as.data.frame(get_pca_ind(pca_res)$coord)
  # # Add clusters obtained using the K-means algorithm
  # ind.coord$cluster <- factor(res.km$cluster)
  # # Add Species groups from the original data sett
  # ind.coord$brainAreas <- windowPCASt$brainAreas
  # 
  
  eigenvalue <- round(get_eigenvalue(pca_res), 1)
  variance.percent <- eigenvalue$variance.percent
  
  # plot<-ggscatter(
  #   ind.coord, x = "Dim.1", y = "Dim.2", 
  #   color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  #   shape = "brainAreas", size = 1.5,  legend = "right", ggtheme = theme_bw()
  # ) + xlab(correspStage[[i]])+ stat_mean(aes(color = cluster), size = 0.5)
  # 
  # plotsKmeans[[i]]<-plot
  
  numPCA<-length(ev[ev>=1]) # Number of selected PCA by Kaiser’s rule
  #numPCA<-ncol(tsne_res) #In case tnse and not PCA
  # 
  # pca_res.ve <- ev/sum(ev)
  # vardf<-data.frame(var=pca_res.ve, PCA=1:length(pca_res.ve))
  # varPCA[[i]]<- ggplot(vardf, aes(PCA, var)) +labs(x="Principal Component",ylab="Proportion of Variance Explained")+
  #   geom_line() +
  #   geom_point()
  
  structDist<-hash()
  for (structure in unique(PCi$topStructure)){
    dfStructure<-PCi %>% filter(topStructure==structure)
    dfOtherStructures<-PCi %>% filter(topStructure!=structure)
    distancesStruct<-c()
    sdStruct<-c()
    for (row in 1:nrow(dfStructure)) {
      pcasStruct<-dfStructure[row,][,1:numPCA]
      # xstruct<-dfStructure[row,]$PC1
      # ystruct<-dfStructure[row,]$PC2
      # zstruct<-dfStructure[row,]$PC3
      nrowsOthers<-nrow(dfOtherStructures)
      listDist<-c()
      for (rowOth in 1:nrowsOthers){
        pcasOtherStruct<-dfOtherStructures[rowOth,][,1:numPCA]
        # xOtherstruct<-dfOtherStructures[rowOth,]$PC1
        # yOtherstruct<-dfOtherStructures[rowOth,]$PC2
        # zOtherstruct<-dfOtherStructures[rowOth,]$PC3
        distance<-dist(matrix(c(as.vector(pcasStruct),as.vector(pcasOtherStruct)),nrow=2,ncol=numPCA),diag=TRUE)
        # distance<-dist(matrix(c(xstruct,ystruct,zstruct,xOtherstruct,yOtherstruct,zOtherstruct),nrow=2,ncol=3),diag=TRUE)
        distance
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

#Testing distances with willcox test

#Initialize dataframe
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


correspStage<-list()
correspStage[[1]]<-"Prenatal"
correspStage[[2]]<-"Infant"
correspStage[[3]]<-"Child"
correspStage[[4]]<-"Adolescent"
correspStage[[5]]<-"Adult"
for (i in 1:5){
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

#wilcoxAkeyDist / wilcoxAkeyPeyDist
#write.xlsx(pairwiseWilcox,"wilcoxAkeyPeyDist.xlsx",col.names=TRUE,row.names=TRUE)

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
#write.csv(wilcoxTests,"wilcoxABAAkeyDist.csv")

willCoxPvals<-as.data.frame(list("brainRegion","window","pvalAVG","log2pvalAVG"))
colnames(willCoxPvals)<-c("brainRegion","window","pvalAVG","log2pvalAVG")
willCoxPvals$pvalAVG<-as.numeric(willCoxPvals$pvalAVG)
willCoxPvals<-willCoxPvals[-1,]
for (col in colnames(willcoxPvaluesAVG)){
  #change the window number
  willCoxPvals<-willCoxPvals %>% add_row(brainRegion=strsplit(col," _ ")[[1]][1],window=strsplit(col," _ ")[[1]][2],pvalAVG=willcoxPvaluesAVG[col][[1]])
}
willCoxPvals$log2pvalAVG<-log2(willCoxPvals$pvalAVG)
colnames(willCoxPvals)<-c("Structure","window","pvalAVG","log2pvalAVG")
plotWill<-ggplot(willCoxPvals, aes(x=window, y=log2pvalAVG, group=Structure)) +
  geom_line(aes(color=Structure),alpha=0.5)+
  geom_point(aes(color=Structure),alpha=0.5)+
  geom_hline(yintercept = log2(0.01), colour="black", size=1.25, alpha=0.5)+
  # scale_x_discrete(breaks = c(2,3,4,5,6,7,8,9),labels=as.character(c("fetal1","fetal2","fetal3","Birth/Infant","Infant/Child","Child","Adolescent","Adult"))) + 
  scale_x_discrete(breaks = c(1,2,3,4,5),labels=as.character(c("Prenatal","Infant","Child","Adolescent","Adult"))) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5))

#Boxplots + log2
library(ggplot2)
library(reshape2)
library(reshape2)
correspStage<-list()
correspStage[[1]]<-"Prenatal"
correspStage[[2]]<-"Infant"
correspStage[[3]]<-"Child"
correspStage[[4]]<-"Adolescent"
correspStage[[5]]<-"Adult"

boxplotsDist<-list()
valoresStruct <- list() 
for (i in 1:5){
  structDist<-structDistWind[[i]]
  for (structure in keys(structDist)){
    valoresStruct[[structure]]<-values(structDist[structure])
  }
  #In a melting pot
  valoresStructdf <- melt(valoresStruct)
  valoresStructdf$Var1<-NULL
  valoresStructdf$L1<-NULL
  
  # prepare a special xlab with the number of obs for each group
  my_xlab <- paste(levels(valoresStructdf$Var1),"\n(N=",table(valoresStructdf$Var1),")",sep="")
  colnames(valoresStructdf) <- c("Structures", "Distance")
  # plot
  boxplotsDist[[i]]<-ggplot(valoresStructdf, aes(x=Structures, y=Distance, fill=Structures)) + geom_boxplot(varwidth = TRUE, alpha=0.5) +
    theme(legend.position="none",axis.text.x = element_blank()) + xlab(correspStage[[i]])
}
supplfig_x<-ggarrange(boxplotsDist[[1]], boxplotsDist[[2]],boxplotsDist[[3]],
          boxplotsDist[[4]],boxplotsDist[[5]],plotWill,
          common.legend = TRUE, legend = "right")
write.csv(willCoxPvals, file="~/AA_investigation/Introgression_deserts/code/new_july/tnse_ABA_distances/wilcox_PeyNotAkeyABA.csv")
# ggsave(file="ABA_Boxplots_DistancesAkeyPey.pdf", width = 11.69, height = 8.27)
ggsave(supplfig_x, file="~/AA_investigation/Introgression_deserts/code/new_july/tnse_ABA_distances/filtered_PeyNotAkeyABAwind_supplfig.pdf", width = 11.69, height = 8.27, units = "in")

# supplfig_x<-ggarrange(plotsKmeans[[1]], plotsKmeans[[2]],plotsKmeans[[3]],
#                       plotsKmeans[[4]],plotsKmeans[[5]],
#                       common.legend = TRUE, legend = "right")
# ggsave(supplfig_x, file="~/AA_investigation/Introgression_deserts/code/new_july/Kmeans_ABA/filtered_PeyNotAkey_supplfig.pdf", width = 11.69, height = 8.27, units = "in")


supplfig_x<-ggarrange(varPCA[[1]], varPCA[[2]],varPCA[[3]],
                      varPCA[[4]],varPCA[[5]],
                      common.legend = TRUE, legend = "right")
ggsave(supplfig_x, file="~/AA_investigation/Introgression_deserts/code/new_july/varPCA/variancePCA_ABA_whole_supplfig.pdf", width = 11.69, height = 8.27, units = "in")


supplfig_x<-ggarrange(plotsPCA[[1]], plotsPCA[[2]],plotsPCA[[3]],
                      plotsPCA[[4]],plotsPCA[[5]],
                      common.legend = TRUE, legend = "right")
ggsave(supplfig_x, file="~/AA_investigation/Introgression_deserts/code/new_july/plotsPCA/PCA_ABA_whole_supplfig.pdf", width = 11.69, height = 8.27, units = "in")

supplfig_x<-ggarrange(jackSplot[[1]], jackSplot[[2]],jackSplot[[3]],
                      jackSplot[[4]],jackSplot[[5]])
ggsave(supplfig_x, file="~/AA_investigation/Introgression_deserts/code/new_july/plotsPCA/JackStraw_ABA_whole_supplfig.pdf", width = 11.69, height = 8.27, units = "in")
