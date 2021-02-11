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

