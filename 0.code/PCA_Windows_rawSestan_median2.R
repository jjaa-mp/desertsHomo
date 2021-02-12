library(matrixStats)
#raw Sestan data is filtered by median > 2 for each gene across all structures/stages
rawmedian2 <- finalrawSestan

rawmedian2[1] <- NULL

#Preparing labels
brainTopStruct<-read.csv("~/raul_tesina/0.code/brainRegionCorresp.csv",sep=":")
names(brainTopStruct) <- c("Regioncode", "FullReg")

#Selecting stage == 2
rawmedian2_w2 <-rawmedian2 %>% filter(Window==2)

win2=merge(brainTopStruct,rawmedian2_w2, by="Regioncode")
pca_resw2 <- prcomp(win2[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw2$x,BrainRegion=win2$Regioncode,topStructure=win2$FullReg)
PCi$topStructure <- str_replace(PCi$topStructure, ("Primary Auditory Cortex|Inferior Temporal Cortex|Orbital Prefrontal Cortex|Superior Temporal Cortex|Primary Visual Cortex|Dorsolateral Prefrontal Cortex|Posterior Inferior Parietal Cortex|Primary Motor Cortex|Medial Prefrontal Cortex|Primary Somatosensory Cortex|Ventrolateral Prefrontal Cortex"), "Neocortex")

w2<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#Selecting stage == 3
rawmedian2_w3 <-rawmedian2 %>% filter(Window==3)

win3=merge(brainTopStruct,rawmedian2_w3, by="Regioncode")
pca_resw3 <- prcomp(win3[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw3$x,BrainRegion=win3$Regioncode,topStructure=win3$FullReg)

w3<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#Selecting stage == 4
rawmedian2_w4 <-rawmedian2 %>% filter(Window==4)

win4=merge(brainTopStruct,rawmedian2_w4, by="Regioncode")
pca_resw4 <- prcomp(win4[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw4$x,BrainRegion=win4$Regioncode,topStructure=win4$FullReg)

w4<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#Selecting stage == 5
rawmedian2_w5 <-rawmedian2 %>% filter(Window==5)

win5=merge(brainTopStruct,rawmedian2_w5, by="Regioncode")
pca_resw5 <- prcomp(win5[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw5$x,BrainRegion=win5$Regioncode,topStructure=win5$FullReg)

w5<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#Selecting stage == 6
rawmedian2_w6 <-rawmedian2 %>% filter(Window==6)

win6=merge(brainTopStruct,rawmedian2_w6, by="Regioncode")
pca_resw6 <- prcomp(win6[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw6$x,BrainRegion=win6$Regioncode,topStructure=win6$FullReg)

w6<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#Selecting stage == 7
rawmedian2_w7 <-rawmedian2 %>% filter(Window==7)

win7=merge(brainTopStruct,rawmedian2_w7, by="Regioncode")
pca_resw7 <- prcomp(win7[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw7$x,BrainRegion=win7$Regioncode,topStructure=win7$FullReg)

w7<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#Selecting stage == 8
rawmedian2_w8 <-rawmedian2 %>% filter(Window==8)

win8=merge(brainTopStruct,rawmedian2_w8, by="Regioncode")
pca_resw8 <- prcomp(win8[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw8$x,BrainRegion=win8$Regioncode,topStructure=win8$FullReg)

w8<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#Selecting stage == 9
rawmedian2_w9 <-rawmedian2 %>% filter(Window==9)

win9=merge(brainTopStruct,rawmedian2_w9, by="Regioncode")
pca_resw9 <- prcomp(win9[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw9$x,BrainRegion=win9$Regioncode,topStructure=win9$FullReg)

w9<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()
ggarrange(w2, w3,w4,w5,w6,w7,w8,w9, labels = c("Fetal1", "Fetal2","Fetal3","Birth/Infan","Infan/Child","Child","Adolescent","Adult"),common.legend = TRUE, legend = "right")
#ggarrange(w2, w3,w4,w5,w6,w7,w8,w9, labels = c("Fetal1", "Fetal2","Fetal3","Birth/Infan","Infan/Child","Child","Adolescent","Adult"),common.legend = TRUE, legend = "right") %>% ggexport(filename = "~/raul_tesina/2.plots/Sestan_raw_median2_PCA/Sestan_raw_media_PCA_stages.pdf")


#Euclidean distance between brain regions. Currently needs to run it for each window
##RawSestan
library(hash)

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
dfDistancesAkeyPey<-data.frame(brainRegion=character(),window=integer(),dist=numeric(),sd=numeric())
structDistw<-as.data.frame(values(structDist))
for (col in colnames(structDistw)){
  #change the window number
  dfDistancesAkeyPey<-dfDistancesAkeyPey %>% add_row(brainRegion=col,window=9,dist=structDistw[col][[1]][1],sd=structDistw[col][[1]][2])
}
write.csv(dfDistancesAkeyPey,"distances_PCAakeyPey9.csv", row.names = FALSE)

file_names <- dir("~/tmp_juan/PC_distances")
final_merged <- do.call(rbind,lapply(file_names,read.csv))


final_mergedTemp<-final_merged
final_mergedTemp$window<-as.character(final_merged$window)

final_mergedTemp$window[final_mergedTemp$window=="2"]<-"2fetal1"
final_mergedTemp$window[final_mergedTemp$window=="3"]<-"3fetal2"
final_mergedTemp$window[final_mergedTemp$window=="4"]<-"4fetal3"
final_mergedTemp$window[final_mergedTemp$window=="5"]<-"5Birth/Infan"
final_mergedTemp$window[final_mergedTemp$window=="6"]<-"6Infan/Child"
final_mergedTemp$window[final_mergedTemp$window=="7"]<-"7Child"
final_mergedTemp$window[final_mergedTemp$window=="8"]<-"8Adolescent"
final_mergedTemp$window[final_mergedTemp$window=="9"]<-"9Adult"

#Plot
library(plotly)

dfDistancesTemp <- final_mergedTemp
fig<-dfDistancesTemp %>% filter(window=="2fetal1") %>%plot_ly(
  x = ~brainRegion, 
  y = ~dist, 
  color=c( 'blue', 'green','grey', 'pink', 'purple', 'red'),
  type = 'bar',
  legendgroup = "A", showlegend=F
) %>% add_trace(
  x = ~brainRegion,y = ~sd,color=c("orange"), showlegend=F) %>%
  layout(xaxis=list(title="fetal2"),yaxis = list(title = 'Dist/sd'), barmode = 'overlay')

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

