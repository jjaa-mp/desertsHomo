library(matrixStats)
#Subsetting raw Sestan data by median > 2 for each gene across all structures/stages
rawmedian2 <- finalrawSestan
rownames(rawmedian2)<- do.call(paste,c(rawmedian2[c("Braincode","Regioncode", "Window")],sep="_"))
rawmedian2[(1:3)] <- NULL
rawmedian2<-rawmedian2 %>% select(which(colMedians(as.matrix(rawmedian2))>2))

rawmedian2 <- rownames_to_column(rawmedian2) 
rawmedian2 <- rawmedian2 %>% 
  separate(rowname, c("Braincode","Regioncode", "Window"))

rawmedian2[1] <- NULL

#Preparing labels
brainTopStruct<-read.csv("~/raul_tesina/0.code/brainRegionCorresp.csv",sep=":")
names(brainTopStruct) <- c("Regioncode", "FullReg")

#Selecting stage == 2
rawmedian2_w2 <-rawmedian2 %>% filter(Window==2)

win2=merge(brainTopStruct,rawmedian2_w2, by="Regioncode")
pca_resw2 <- prcomp(win2[-(1:3)], scale. = TRUE)

PCi<-data.frame(pca_resw2$x,BrainRegion=win2$Regioncode,topStructure=win2$FullReg)

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
