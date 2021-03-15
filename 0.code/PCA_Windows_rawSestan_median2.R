#rawSestan
```{r}
load("~/finalraw2_inputPC.RData")
rawmedian2 <- finalraw2 #Input from jm_Allen or load("finalraw2_inputPC.RData")

rawmedian2[1] <- NULL

structDistWind<-list()

for (i in 2:9){
  
  windowPCA<- rawmedian2 %>% filter(Window==i)

  windowPCA$Window<-NULL
  
  brainTopStruct<-read.csv("~/raul_tesina/0.code/brainRegionCorresp.csv",sep=":")
  colnames(brainTopStruct) <- c("Regioncode","FullReg")

  windowPCASt=merge(brainTopStruct,windowPCA,by="Regioncode")
  
  pca_res <- prcomp(windowPCASt[,-1][,-1], scale. = TRUE)
  
  # PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$topStructure) 
  PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$FullReg)
  PCi$topStructure <- str_replace(PCi$topStructure, ("Primary Auditory Cortex|Inferior Temporal Cortex|Orbital Prefrontal Cortex|Superior Temporal Cortex|Primary Visual Cortex|Dorsolateral Prefrontal Cortex|Posterior Inferior Parietal Cortex|Primary Motor Cortex|Medial Prefrontal Cortex|Primary Somatosensory Cortex|Ventrolateral Prefrontal Cortex"), "Neocortex")
 
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

#Boxplots
correspStage<-list()
correspStage[[2]]<-"Fetal_1"
correspStage[[3]]<-"Fetal_2"
correspStage[[4]]<-"Fetal_3"
correspStage[[5]]<-"Birth/Inf"
correspStage[[6]]<-"Infan_Child"
correspStage[[7]]<-"Child"
correspStage[[8]]<-"Adolescent"
correspStage[[9]]<-"Adult"

boxplotsDist<-list()
valoresStruct <- list() 
for (i in 2:9){
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

#save(boxplotsDist, file = "boxplotsFilteredRawSestan.R")

#BOXPLOTS:
#load("~/tmp_wilcoxtests/boxplotsFilteredRawSestan.R")

boxplotsDist[[2]]
ggarrange(boxplotsDist[[2]], boxplotsDist[[3]],boxplotsDist[[4]],
          boxplotsDist[[5]],boxplotsDist[[6]],boxplotsDist[[7]],
          boxplotsDist[[8]],boxplotsDist[[9]],
          common.legend = TRUE, legend = "right")

#Pairwise wilcox tests
wilcoxTests<-data.frame(matrix(ncol=48,nrow = 6))
wilcoxTestsCol<-c()
for (wind in 2:9){
  for (str in keys(structDist)){
    wilcoxTestsCol<-append(wilcoxTestsCol,paste(str,"_",wind))
  }
}
colnames(wilcoxTests)<-wilcoxTestsCol
rownames(wilcoxTests)<-keys(structDist)
pairwiseWilcox<-list()
for (i in 2:9){
  structDist<-structDistWind[[i]]
  structs<-c()
  valuesDist<-c()
  for(str in keys(structDist)){
    structs<-append(structs,rep(str,length(values(structDist[str]))))
    valuesDist<-append(valuesDist,values(structDist[str]))
  }
  pairwiseWilcox[[correspStage[[i]]]]<-as.data.frame(pairwise.wilcox.test(valuesDist,structs, p.adj = "bonf")$p.value)
}
#save(pairwiseWilcox, "pairwiseWilcox_filteredRawSestan.RData")
#write.xlsx(pairwiseWilcox,"pairwiseWilcoxFilteredRawSestan.xlsx",col.names=TRUE,row.names=TRUE)

#load("~/tmp_wilcoxtests/pairwiseWilcox_filteredRawSestan.RData")

colMeans(wilcoxTests, na.rm = TRUE)

wilcoxTests<-data.frame(matrix(ncol=48,nrow = 6))
wilcoxTestsCol<-c()
for (wind in 2:9){
  for (str in keys(structDist)){
    wilcoxTestsCol<-append(wilcoxTestsCol,paste(str,"_",wind))
  }
}
colnames(wilcoxTests)<-wilcoxTestsCol
rownames(wilcoxTests)<-keys(structDist)
for (i in 2:9){
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
#write.csv(wilcoxTests,"wilcoxAkeyDist.csv")

willCoxPvals<-as.data.frame(list("Structure","window","log2pvalAVG"))
colnames(willCoxPvals)<-c("Structure","window","log2pvalAVG")
willCoxPvals$log2pvalAVG<-as.numeric(willCoxPvals$log2pvalAVG)
willCoxPvals<-willCoxPvals[-1,]
for (col in colnames(willcoxPvaluesAVG)){
  #change the window number
  willCoxPvals<-willCoxPvals %>% add_row(Structure=strsplit(col," _ ")[[1]][1],window=strsplit(col," _ ")[[1]][2],log2pvalAVG=willcoxPvaluesAVG[col][[1]])
}
willCoxPvals$log2pvalAVG<-log2(willCoxPvals$log2pvalAVG)


willCoxPvals$window <- str_replace(willCoxPvals$window, "2", "Fetal_1")
willCoxPvals$window <- str_replace(willCoxPvals$window, "3", "Fetal_2")
willCoxPvals$window <- str_replace(willCoxPvals$window, "4", "Fetal_3")
willCoxPvals$window <- str_replace(willCoxPvals$window, "5", "Birth/Inf")
willCoxPvals$window <- str_replace(willCoxPvals$window, "6", "Inf/Child")
willCoxPvals$window <- str_replace(willCoxPvals$window, "7", "Child")
willCoxPvals$window <- str_replace(willCoxPvals$window, "8", "Adolescence")
willCoxPvals$window <- str_replace(willCoxPvals$window, "9", "Adult")

willCoxPvals$window = factor(willCoxPvals$window, 
    levels=c("Fetal_1","Fetal_2","Fetal_3","Birth/Inf","Inf/Child", "Child","Adolescence", "Adult"))

#save(willCoxPvals, file="ggplot_wilcox_filteredRawSestan")
#load("ggplot_wilcox_filteredRawSestan")

#write.csv(willCoxPvals, file="~/raul_tesina/1.data/distances_pvalues/wilcox_filteredRawSestan.csv")

##PREPARING SUPPL FIGURE
#load("~/tmp_wilcoxtests/ggplot_wilcox_filteredRawSestan")
subfigx <- ggplot(willCoxPvals, aes(x=window, y=log2pvalAVG, group=Structure)) +
  geom_line(aes(color=Structure))+
  geom_point(aes(color=Structure))+xlab("")+ylab("log2pvalAVG")+theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 45,  hjust = 1), legend.position = "right")+ggtitle("Pairwise Wilcox test")+geom_hline(yintercept = log2(0.01), colour="black", size=1.25, alpha=0.5)

#load("~/tmp_wilcoxtests/boxplotsFilteredRawSestan.R")
ggarrange(boxplotsDist[[2]], boxplotsDist[[3]],boxplotsDist[[4]],
          boxplotsDist[[5]],boxplotsDist[[6]],boxplotsDist[[7]],
          boxplotsDist[[8]],boxplotsDist[[9]],
          common.legend = TRUE, legend = "right")

supplfig_x <- ggarrange(boxplotsDist[[2]], boxplotsDist[[3]],boxplotsDist[[4]],
          boxplotsDist[[5]],boxplotsDist[[6]],boxplotsDist[[7]],
          boxplotsDist[[8]],boxplotsDist[[9]], subfigx,
          common.legend = TRUE, legend = "right")

#ggsave(supplfig_x, file="~/raul_tesina/2.plots/Sestan_raw_filtered_median2_trajectories/filtered_rawSestan_supplfig.pdf", width = 11.69, height = 8.27, units = "in")
```

#AKEY
```{r}
akeySestan <- finalakeySestan1 #Input from jm_Allen

akeySestan[1] <- NULL

structDistWind<-list()

for (i in 2:9){
  
#akeySestan <- from jm_Allen finalraw2 inputPC
  windowPCA<- akeySestan %>% filter(Window==i)

  windowPCA$Window<-NULL
  
  brainTopStruct<-read.csv("~/raul_tesina/0.code/brainRegionCorresp.csv",sep=":")
  colnames(brainTopStruct) <- c("Regioncode","FullReg")

  windowPCASt=merge(brainTopStruct,windowPCA,by="Regioncode")
  
  pca_res <- prcomp(windowPCASt[,-1][,-1], scale. = TRUE)
  
  # PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$topStructure) 
  PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$FullReg)
  PCi$topStructure <- str_replace(PCi$topStructure, ("Primary Auditory Cortex|Inferior Temporal Cortex|Orbital Prefrontal Cortex|Superior Temporal Cortex|Primary Visual Cortex|Dorsolateral Prefrontal Cortex|Posterior Inferior Parietal Cortex|Primary Motor Cortex|Medial Prefrontal Cortex|Primary Somatosensory Cortex|Ventrolateral Prefrontal Cortex"), "Neocortex")
 
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

#Boxplots
correspStage<-list()
correspStage[[2]]<-"Fetal_1"
correspStage[[3]]<-"Fetal_2"
correspStage[[4]]<-"Fetal_3"
correspStage[[5]]<-"Birth/Inf"
correspStage[[6]]<-"Infan_Child"
correspStage[[7]]<-"Child"
correspStage[[8]]<-"Adolescent"
correspStage[[9]]<-"Adult"

boxplotsDist<-list()
valoresStruct <- list() 
for (i in 2:9){
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

#save(boxplotsDist, file = "boxplotsFilteredRawSestan.R")

#BOXPLOTS:
#load("~/tmp_wilcoxtests/boxplotsFilteredRawSestan.R")

boxplotsDist[[2]]
ggarrange(boxplotsDist[[2]], boxplotsDist[[3]],boxplotsDist[[4]],
          boxplotsDist[[5]],boxplotsDist[[6]],boxplotsDist[[7]],
          boxplotsDist[[8]],boxplotsDist[[9]],
          common.legend = TRUE, legend = "right")

#Pairwise wilcox tests
wilcoxTests<-data.frame(matrix(ncol=48,nrow = 6))
wilcoxTestsCol<-c()
for (wind in 2:9){
  for (str in keys(structDist)){
    wilcoxTestsCol<-append(wilcoxTestsCol,paste(str,"_",wind))
  }
}
colnames(wilcoxTests)<-wilcoxTestsCol
rownames(wilcoxTests)<-keys(structDist)
pairwiseWilcox<-list()
for (i in 2:9){
  structDist<-structDistWind[[i]]
  structs<-c()
  valuesDist<-c()
  for(str in keys(structDist)){
    structs<-append(structs,rep(str,length(values(structDist[str]))))
    valuesDist<-append(valuesDist,values(structDist[str]))
  }
  pairwiseWilcox[[correspStage[[i]]]]<-as.data.frame(pairwise.wilcox.test(valuesDist,structs, p.adj = "bonf")$p.value)
}
#save(pairwiseWilcox, "pairwiseWilcox_filteredRawSestan.RData")
#write.xlsx(pairwiseWilcox,"pairwiseWilcoxFilteredRawSestan.xlsx",col.names=TRUE,row.names=TRUE)

#load("~/tmp_wilcoxtests/pairwiseWilcox_filteredRawSestan.RData")

colMeans(wilcoxTests, na.rm = TRUE)

wilcoxTests<-data.frame(matrix(ncol=48,nrow = 6))
wilcoxTestsCol<-c()
for (wind in 2:9){
  for (str in keys(structDist)){
    wilcoxTestsCol<-append(wilcoxTestsCol,paste(str,"_",wind))
  }
}
colnames(wilcoxTests)<-wilcoxTestsCol
rownames(wilcoxTests)<-keys(structDist)
for (i in 2:9){
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
#write.csv(wilcoxTests,"wilcoxAkeyDist.csv")

willCoxPvals<-as.data.frame(list("Structure","window","log2pvalAVG"))
colnames(willCoxPvals)<-c("Structure","window","log2pvalAVG")
willCoxPvals$log2pvalAVG<-as.numeric(willCoxPvals$log2pvalAVG)
willCoxPvals<-willCoxPvals[-1,]
for (col in colnames(willcoxPvaluesAVG)){
  #change the window number
  willCoxPvals<-willCoxPvals %>% add_row(Structure=strsplit(col," _ ")[[1]][1],window=strsplit(col," _ ")[[1]][2],log2pvalAVG=willcoxPvaluesAVG[col][[1]])
}
willCoxPvals$log2pvalAVG<-log2(willCoxPvals$log2pvalAVG)


willCoxPvals$window <- str_replace(willCoxPvals$window, "2", "Fetal_1")
willCoxPvals$window <- str_replace(willCoxPvals$window, "3", "Fetal_2")
willCoxPvals$window <- str_replace(willCoxPvals$window, "4", "Fetal_3")
willCoxPvals$window <- str_replace(willCoxPvals$window, "5", "Birth/Inf")
willCoxPvals$window <- str_replace(willCoxPvals$window, "6", "Inf/Child")
willCoxPvals$window <- str_replace(willCoxPvals$window, "7", "Child")
willCoxPvals$window <- str_replace(willCoxPvals$window, "8", "Adolescence")
willCoxPvals$window <- str_replace(willCoxPvals$window, "9", "Adult")

willCoxPvals$window = factor(willCoxPvals$window, 
    levels=c("Fetal_1","Fetal_2","Fetal_3","Birth/Inf","Inf/Child", "Child","Adolescence", "Adult"))

ggplot(willCoxPvals, aes(x=window, y=log2pvalAVG, group=Structure)) +
  geom_line(aes(color=Structure))+
  geom_point(aes(color=Structure))+xlab("")+theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 45,  hjust = 1))
pvalAKEY <- ggplot(willCoxPvals, aes(x=window, y=log2pvalAVG, group=Structure)) +
  geom_line(aes(color=Structure))+
  geom_point(aes(color=Structure))+xlab("")+theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 45,  hjust = 1))+geom_hline(yintercept = log2(0.01), colour="black", size=1.25, alpha=0.5)
#write.csv(willCoxPvals, file="~/raul_tesina/1.data/distances_pvalues/wilcox_akeySestan.csv")
```

#AKEPEY
```{r}
akeypeySestan <- finalakeypeySestan #Input from jm_Allen or load("finalraw2_inputPC.RData")

akeypeySestan[1] <- NULL

structDistWind<-list()

for (i in 2:9){
  windowPCA<- akeypeySestan %>% filter(Window==i)

  windowPCA$Window<-NULL
  
  brainTopStruct<-read.csv("~/raul_tesina/0.code/brainRegionCorresp.csv",sep=":")
  colnames(brainTopStruct) <- c("Regioncode","FullReg")

  windowPCASt=merge(brainTopStruct,windowPCA,by="Regioncode")
  
  pca_res <- prcomp(windowPCASt[,-1][,-1], scale. = TRUE)
  
  # PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$topStructure) 
  PCi<-data.frame(pca_res$x,BrainRegion=windowPCASt$Regioncode,topStructure=windowPCASt$FullReg)
  PCi$topStructure <- str_replace(PCi$topStructure, ("Primary Auditory Cortex|Inferior Temporal Cortex|Orbital Prefrontal Cortex|Superior Temporal Cortex|Primary Visual Cortex|Dorsolateral Prefrontal Cortex|Posterior Inferior Parietal Cortex|Primary Motor Cortex|Medial Prefrontal Cortex|Primary Somatosensory Cortex|Ventrolateral Prefrontal Cortex"), "Neocortex")
 
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

#Boxplots
correspStage<-list()
correspStage[[2]]<-"Fetal_1"
correspStage[[3]]<-"Fetal_2"
correspStage[[4]]<-"Fetal_3"
correspStage[[5]]<-"Birth/Inf"
correspStage[[6]]<-"Infan_Child"
correspStage[[7]]<-"Child"
correspStage[[8]]<-"Adolescent"
correspStage[[9]]<-"Adult"

boxplotDistakeypey<-list()
valoresStruct <- list() 
for (i in 2:9){
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
  boxplotDistakeypey[[i]]<-ggplot(valoresStructdf, aes(x=Structures, y=Distance, fill=Structures)) + geom_boxplot(varwidth = TRUE, alpha=0.5) +
  theme(legend.position="none",axis.text.x = element_blank()) + xlab(correspStage[[i]])
}

#save(boxplotDistakeypey, file = "boxplotsFilteredRawSestan.R")

#BOXPLOTS:
#load("~/tmp_wilcoxtests/boxplotsFilteredRawSestan.R")

boxplotDistakeypey[[2]]
ggarrange(boxplotDistakeypey[[2]], boxplotDistakeypey[[3]],boxplotDistakeypey[[4]],
          boxplotDistakeypey[[5]],boxplotDistakeypey[[6]],boxplotDistakeypey[[7]],
          boxplotDistakeypey[[8]],boxplotDistakeypey[[9]],
          common.legend = TRUE, legend = "right")

#Pairwise wilcox tests
wilcoxTests<-data.frame(matrix(ncol=48,nrow = 6))
wilcoxTestsCol<-c()
for (wind in 2:9){
  for (str in keys(structDist)){
    wilcoxTestsCol<-append(wilcoxTestsCol,paste(str,"_",wind))
  }
}
colnames(wilcoxTests)<-wilcoxTestsCol
rownames(wilcoxTests)<-keys(structDist)
pairwiseWilcox<-list()
for (i in 2:9){
  structDist<-structDistWind[[i]]
  structs<-c()
  valuesDist<-c()
  for(str in keys(structDist)){
    structs<-append(structs,rep(str,length(values(structDist[str]))))
    valuesDist<-append(valuesDist,values(structDist[str]))
  }
  pairwiseWilcox[[correspStage[[i]]]]<-as.data.frame(pairwise.wilcox.test(valuesDist,structs, p.adj = "bonf")$p.value)
}
#save(pairwiseWilcox, "pairwiseWilcox_filteredRawSestan.RData")
#write.xlsx(pairwiseWilcox,"pairwiseWilcoxFilteredRawSestan.xlsx",col.names=TRUE,row.names=TRUE)

#load("~/tmp_wilcoxtests/pairwiseWilcox_filteredRawSestan.RData")

colMeans(wilcoxTests, na.rm = TRUE)

wilcoxTests<-data.frame(matrix(ncol=48,nrow = 6))
wilcoxTestsCol<-c()
for (wind in 2:9){
  for (str in keys(structDist)){
    wilcoxTestsCol<-append(wilcoxTestsCol,paste(str,"_",wind))
  }
}
colnames(wilcoxTests)<-wilcoxTestsCol
rownames(wilcoxTests)<-keys(structDist)
for (i in 2:9){
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
#write.csv(wilcoxTests,"wilcoxAkeyDist.csv")

willCoxPvals<-as.data.frame(list("Structure","window","log2pvalAVG"))
colnames(willCoxPvals)<-c("Structure","window","log2pvalAVG")
willCoxPvals$log2pvalAVG<-as.numeric(willCoxPvals$log2pvalAVG)
willCoxPvals<-willCoxPvals[-1,]
for (col in colnames(willcoxPvaluesAVG)){
  #change the window number
  willCoxPvals<-willCoxPvals %>% add_row(Structure=strsplit(col," _ ")[[1]][1],window=strsplit(col," _ ")[[1]][2],log2pvalAVG=willcoxPvaluesAVG[col][[1]])
}
willCoxPvals$log2pvalAVG<-log2(willCoxPvals$log2pvalAVG)


willCoxPvals$window <- str_replace(willCoxPvals$window, "2", "Fetal_1")
willCoxPvals$window <- str_replace(willCoxPvals$window, "3", "Fetal_2")
willCoxPvals$window <- str_replace(willCoxPvals$window, "4", "Fetal_3")
willCoxPvals$window <- str_replace(willCoxPvals$window, "5", "Birth/Inf")
willCoxPvals$window <- str_replace(willCoxPvals$window, "6", "Inf/Child")
willCoxPvals$window <- str_replace(willCoxPvals$window, "7", "Child")
willCoxPvals$window <- str_replace(willCoxPvals$window, "8", "Adolescence")
willCoxPvals$window <- str_replace(willCoxPvals$window, "9", "Adult")

willCoxPvals$window = factor(willCoxPvals$window, 
    levels=c("Fetal_1","Fetal_2","Fetal_3","Birth/Inf","Inf/Child", "Child","Adolescence", "Adult"))


ggplot(willCoxPvals, aes(x=window, y=log2pvalAVG, group=Structure)) +
  geom_line(aes(color=Structure))+
  geom_point(aes(color=Structure))+xlab("")+theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 45,  hjust = 1))
pvalAKEYPEY <- ggplot(willCoxPvals, aes(x=window, y=log2pvalAVG, group=Structure)) +
  geom_line(aes(color=Structure))+
  geom_point(aes(color=Structure))+xlab("")+theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 45,  hjust = 1))+geom_hline(yintercept = log2(0.01), colour="black", size=1.25, alpha=0.5)
#write.csv(willCoxPvals, file="~/raul_tesina/1.data/distances_pvalues/wilcox_akeypeySestan.csv")
```

#Table p-values
```{r}
rawp <- read.csv("~/raul_tesina/1.data/distances_pvalues/wilcox_filteredRawSestan.csv")
colnames(rawp) <- c("Structure", "window", "log2(pvalAVG)") 
rawp$pvalAVG <- 2**(rawp$`log2(pvalAVG)`)
order <- c("Structure", "window", "pvalAVG", "log2(pvalAVG)")
rawp <- rawp[, order]

akeyp <- read.csv("~/raul_tesina/1.data/distances_pvalues/wilcox_akeySestan.csv", row.names = 1)
colnames(akeyp) <- c("Structure", "window", "log2(pvalAVG)") 
akeyp$pvalAVG <- 2**(akeyp$`log2(pvalAVG)`)
order <- c("Structure", "window", "pvalAVG", "log2(pvalAVG)")
akeyp <- akeyp[, order]

akeypeyp <- read.csv("~/raul_tesina/1.data/distances_pvalues/wilcox_akeypeySestan.csv", row.names = 1)
colnames(akeypeyp) <- c("Structure", "window", "log2(pvalAVG)") 
akeypeyp$pvalAVG <- 2**(akeypeyp$`log2(pvalAVG)`)
order <- c("Structure", "window", "pvalAVG", "log2(pvalAVG)")
akeypeyp <- akeypeyp[, order]

# write.csv(rawp, file="~/raul_tesina/1.data/distances_pvalues/wilcox_filteredRawSestan.csv")
# write.csv(akeyp, file="~/raul_tesina/1.data/distances_pvalues/wilcox_akeySestan.csv")
# write.csv(akeypeyp, file="~/raul_tesina/1.data/distances_pvalues/wilcox_akeypeySestan.csv")
```

#Figure p-values
```{r}
pvala<-arrangeGrob(pvalAKEY, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top")))
pvalb<-arrangeGrob(pvalAKEYPEY, top = textGrob("B", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top")))
pvals_inter<-grid.arrange(pvala, pvalb)

ggsave(file="~/raul_tesina/2.plots/Sestan_AkeyPey_log2_median/sf_pvals_ak_akpey.pdf", pvals_inter, width = 11.69, height = 8.27, units = "in")
```