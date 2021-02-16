#Euclidean distance between brain regions
library(hash)
#AkeyPey
structDist<-hash()
for (structure in unique(PCi$topStructure)){
  
  dfStructure<-PCi %>% filter(topStructure==structure)
  dfOtherStructures<-PCi %>% filter(topStructure!=structure)
  meansStruct<-c()
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
    meansStruct<-append(meansStruct,mean(listDist))
    sdStruct<-append(sdStruct,sd(listDist))
    print(paste("Distance ",structure,": mean->",mean(listDist),"; sd->",sd(listDist)))
  }
  structDist[[structure]]<-c(mean(meansStruct),mean(sdStruct))
}

#initialize dataframe
#dfDistancesAkeyPey<-data.frame(brainRegion=character(),window=integer(),dist=numeric(),sd=numeric())
structDistw<-as.data.frame(values(structDist))
for (col in colnames(structDistw)){
  #change the window number
  dfDistancesAkeyPey<-dfDistancesAkeyPey %>% add_row(brainRegion=col,window=9,dist=structDistw[col][[1]][1],sd=structDistw[col][[1]][2])
}
dfDistancesAkeyPey
write.csv(dfDistancesAkeyPey,"distances_PCAakeyPey.csv")
dfDistancesAkeyPey
df$Ozone[is.na(df$Ozone<-mean(df$Ozone, na.rm = T)
dfDistancesAkeyPeyTemp<-dfDistancesAkeyPey
dfDistancesAkeyPeyTemp$window<-as.character(dfDistancesAkeyPey$window)

dfDistancesAkeyPeyTemp$window[dfDistancesAkeyPeyTemp$window=="2"]<-"2fetal1"
dfDistancesAkeyPeyTemp$window[dfDistancesAkeyPeyTemp$window=="3"]<-"3fetal2"
dfDistancesAkeyPeyTemp$window[dfDistancesAkeyPeyTemp$window=="4"]<-"4fetal3"
dfDistancesAkeyPeyTemp$window[dfDistancesAkeyPeyTemp$window=="5"]<-"5Birth/Infan"
dfDistancesAkeyPeyTemp$window[dfDistancesAkeyPeyTemp$window=="6"]<-"6Infan/Child"
dfDistancesAkeyPeyTemp$window[dfDistancesAkeyPeyTemp$window=="7"]<-"7Child"
dfDistancesAkeyPeyTemp$window[dfDistancesAkeyPeyTemp$window=="8"]<-"8Adolescent"
dfDistancesAkeyPeyTemp$window[dfDistancesAkeyPeyTemp$window=="9"]<-"9Adult"
#"Birth/Infan","Infan/Child","Child","Adolescent","Adult"


#Akey
structDist<-hash()
for (structure in unique(PCi$topStructure)){
  
  dfStructure<-PCi %>% filter(topStructure==structure)
  dfOtherStructures<-PCi %>% filter(topStructure!=structure)
  meansStruct<-c()
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
    meansStruct<-append(meansStruct,mean(listDist))
    sdStruct<-append(sdStruct,sd(listDist))
    print(paste("Distance ",structure,": mean->",mean(listDist),"; sd->",sd(listDist)))
  }
  structDist[[structure]]<-c(mean(meansStruct),mean(sdStruct))
}

#initialize dataframe
#dfDistancesAkey<-data.frame(brainRegion=character(),window=integer(),dist=numeric(),sd=numeric())
structDistw<-as.data.frame(values(structDist))
for (col in colnames(structDistw)){
  #change the window number
  dfDistancesAkey<-dfDistancesAkey %>% add_row(brainRegion=col,window=9,dist=structDistw[col][[1]][1],sd=structDistw[col][[1]][2])
}
write.csv(dfDistancesAkey,"distances_PCAakey.csv")


#Visualization distances
library(plotly)
dfDistancesTemp<-dfDistancesAkey #dfDistancesAkeyPey/dfDistancesAkey
dfDistancesTemp$window[dfDistancesTemp$window=="2"]<-"2fetal1"
dfDistancesTemp$window[dfDistancesTemp$window=="3"]<-"3fetal2"
dfDistancesTemp$window[dfDistancesTemp$window=="4"]<-"4fetal3"
dfDistancesTemp$window[dfDistancesTemp$window=="5"]<-"5Birth/Infan"
dfDistancesTemp$window[dfDistancesTemp$window=="6"]<-"6Infan/Child"
dfDistancesTemp$window[dfDistancesTemp$window=="7"]<-"7Child"
dfDistancesTemp$window[dfDistancesTemp$window=="8"]<-"8Adolescent"
dfDistancesTemp$window[dfDistancesTemp$window=="9"]<-"9Adult"
fig<-dfDistancesTemp %>% filter(window=="2fetal1") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="fetal1"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')
fig2<-dfDistancesTemp %>% filter(window=="3fetal2") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="fetal2"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')
fig3<-dfDistancesTemp %>% filter(window=="4fetal3") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="fetal3"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')

fig4<-dfDistancesTemp %>% filter(window=="5Birth/Infan") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="Birth/Infan"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')

fig5<-dfDistancesTemp %>% filter(window=="6Infan/Child") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="Infan/Child"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')

fig6<-dfDistancesTemp %>% filter(window=="7Child") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="Child"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')

fig7<-dfDistancesTemp %>% filter(window=="8Adolescent") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="Adolescent"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')

fig8<-dfDistancesTemp %>% filter(window=="9Adult") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="Adult"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')


subplot(fig,fig2,fig3,fig4,fig5,fig6,fig7,fig8,shareX = T) %>% layout(yaxis=list(title = 'Dist/sd'))
              
#Tests distances + Boxplots
#Akey or AkeyPey distances tests

structDistWind<-list()
for (i in 2:9){
  
#CHANGE:logakeypeyPCA / logakeyPCA
  windowPCA<-logakeyPCA %>% filter(Window==i)
  
  windowPCA$Window<-NULL
  
  brainTopStruct<-read.csv("brainRegionCorresp.csv",sep=":")
  windowPCASt=merge(brainTopStruct,windowPCA,by="Regioncode")
  
  pca_res <- prcomp(windowPCASt[,-1][,-1][,-1][,-1], scale. = TRUE)
  
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

#Boxplots
correspStage<-list()
correspStage[[2]]<-"fetal1"
correspStage[[3]]<-"fetal2"
correspStage[[4]]<-"fetal3"
correspStage[[5]]<-"Birth_Infan"
correspStage[[6]]<-"Infan_Child"
correspStage[[7]]<-"Child"
correspStage[[8]]<-"Adolescent"
correspStage[[9]]<-"Adult"
boxplotsDist<-list()
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
colnames(valoresStructdf) <- c("names", "value")
# plot
boxplotsDist[[i]]<-ggplot(valoresStructdf, aes(x=names, y=value, fill=names)) +
  geom_boxplot(varwidth = TRUE, alpha=0.2) +
  theme(legend.position="none",axis.text.x = element_blank()) + xlab(correspStage[[i]])

}
boxplotsDist[[2]]
ggarrange(boxplotsDist[[2]], boxplotsDist[[3]],boxplotsDist[[4]],
          boxplotsDist[[5]],boxplotsDist[[6]],boxplotsDist[[7]],
          boxplotsDist[[8]],boxplotsDist[[9]],
          common.legend = TRUE, legend = "right")
ggsave(file="Sestan_Boxplots_DistancesAkey.pdf", width = 11.69, height = 8.27)
               
#Tests we used Wilcox, other tests can be performed
#Initialize dataframe
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
write.xlsx(pairwiseWilcox,"wilcoxAkeyPeyDist.xlsx",col.names=TRUE,row.names=TRUE)
colMeans(wilcoxTests, na.rm = TRUE)
# wilcoxAkeyDist.csv/wilcoxAkeyPeyDist.csv
write.csv(wilcoxTests,"wilcoxAkeyPeyDist.csv")
               

