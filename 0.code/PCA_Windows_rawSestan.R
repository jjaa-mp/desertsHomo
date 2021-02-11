#Raw Sestan data (without filtering by median > 2)

finalrawSestan ##Input from jm_Allen.R

finalrawSestan$Braincode <- NULL
finalrawSestan$Label<- paste(finalrawSestan$Regioncode, finalrawSestan$Window, sep="_")

rawpc <- finalrawSestan
rawpc$Label <- NULL

#Preparing labels
brainTopStruct<-read.csv("~/raul_tesina/0.code/brainRegionCorresp.csv",sep=":")
names(brainTopStruct) <- c("Regioncode", "FullReg")

#2 
rawpc2 <-rawpc %>% filter(Window==2)

window2=merge(brainTopStruct,rawpc2, by="Regioncode")

pca_res2 <- prcomp(window2[,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res2$x,BrainRegion=window2$Regioncode,topStructure=window2$FullReg)

w2<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#3
rawpc3 <-rawpc %>% filter(Window==3)
rawpc3$Window <- NULL

window3=merge(brainTopStruct,rawpc3, by="Regioncode")

pca_res3 <- prcomp(window3[,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res3$x,BrainRegion=window3$Regioncode,topStructure=window3$FullReg)

w3<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#4
rawpc4 <-rawpc %>% filter(Window==4)
rawpc4$Window <- NULL

window4=merge(brainTopStruct,rawpc4, by="Regioncode")

pca_res4 <- prcomp(window4[,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res4$x,BrainRegion=window4$Regioncode,topStructure=window4$FullReg)

w4<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#5
rawpc5 <-rawpc %>% filter(Window==5)
rawpc5$Window <- NULL

window5=merge(brainTopStruct,rawpc5, by="Regioncode")

pca_res5 <- prcomp(window5[,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res5$x,BrainRegion=window5$Regioncode,topStructure=window5$FullReg)

w5<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()
#6
rawpc6 <-rawpc %>% filter(Window==6)
rawpc6$Window <- NULL

window6=merge(brainTopStruct,rawpc6, by="Regioncode")

pca_res6 <- prcomp(window6[,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res6$x,BrainRegion=window6$Regioncode,topStructure=window6$FullReg)

w6<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#7
rawpc7 <-rawpc %>% filter(Window==7)
rawpc7$Window <- NULL

window7=merge(brainTopStruct,rawpc7, by="Regioncode")

pca_res7 <- prcomp(window7[,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res7$x,BrainRegion=window7$Regioncode,topStructure=window7$FullReg)

w7<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#8
rawpc8 <-rawpc %>% filter(Window==8)
rawpc8$Window <- NULL

window8=merge(brainTopStruct,rawpc8, by="Regioncode")

pca_res8 <- prcomp(window8[,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res8$x,BrainRegion=window8$Regioncode,topStructure=window8$FullReg)

w8<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#9
rawpc9 <-rawpc %>% filter(Window==9)
rawpc9$Window <- NULL

window9=merge(brainTopStruct,rawpc9, by="Regioncode")

pca_res9 <- prcomp(window9[,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res9$x,BrainRegion=window9$Regioncode,topStructure=window9$FullReg)

w9<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  #scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()

#Plotting all stages
ggarrange(w2, w3,w4,w5,w6,w7,w8,w9, labels = c("Fetal1", "Fetal2","Fetal3","Birth/Infan","Infan/Child","Child","Adolescent","Adult"),
          common.legend = TRUE, legend = "right")
