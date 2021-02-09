#PCA AkeyPey

akeyPeyNorm=read.csv("logakeypey.csv")
akeyPeyNorm<-akeyPeyNorm[,-1]
akeyNorm=read.csv("logakey.csv")

logakeypeyPCA<-akeyPeyNorm
logakeypeyPCA$X<-NULL
logakeypeyPCA$Braincode<-NULL

str(logakeypeyPCA)


windowAkeyPeyPCA<-logakeypeyPCA %>% filter(Window==9)
windowAkeyPeyPCA$Window<-NULL

brainTopStruct<-read.csv("brainRegionCorresp.csv",sep=":")
windowAkeyPeyPCASt=merge(brainTopStruct,windowAkeyPeyPCA,by="Regioncode")

pca_res <- prcomp(windowAkeyPeyPCASt[,-1][,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res$x,BrainRegion=windowAkeyPeyPCASt$Regioncode,topStructure=windowAkeyPeyPCASt$topStructure)
w9<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()#+
  #   theme(legend.position = "none")
# "embrionic","fetal1", "fetal2","fetal3","Birth/Infan","Infan/Child","Child","Adolescent","Adult"]
w2
# Extract the legend. Returns a gtable
leg <- get_legend(w2)

# Convert to a ggplot and print
legend<-as_ggplot(leg)
legend
ggarrange(w2, w3,w4,w5,w6,w7,w8,w9, labels = c("Fetal1", "Fetal2","Fetal3","Birth/Infan","Infan/Child","Child","Adolescent","Adult"),
          common.legend = TRUE, legend = "right")

ggsave(file="Sestan_PCA_WindowsAkeyPey.pdf", width = 11.69, height = 8.27)

#PCA Akey
logakeyPCA<-akeyNorm
logakeyPCA$X<-NULL
logakeyPCA$Braincode<-NULL

windowAkeyPCA<-logakeyPCA %>% filter(Window==9)
windowAkeyPCA$Window<-NULL

brainTopStruct<-read.csv("brainRegionCorresp.csv",sep=":")
windowAkeyPCASt=merge(brainTopStruct,windowAkeyPCA,by="Regioncode")

pca_res <- prcomp(windowAkeyPCASt[,-1][,-1][,-1], scale. = TRUE)

PCi<-data.frame(pca_res$x,BrainRegion=windowAkeyPCASt$Regioncode,topStructure=windowAkeyPCASt$topStructure)
w9Akey<-ggplot(PCi,aes(x=PC1,y=PC2,col=topStructure))+
  geom_point(size=3,alpha=1)+ #Size and alpha just for fun
  scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe'))+ #your colors here
  theme_classic()#+
#   theme(legend.position = "none")
# "embrionic","fetal1", "fetal2","fetal3","Birth/Infan","Infan/Child","Child","Adolescent","Adult"]

ggarrange(w2Akey, w3Akey,w4Akey,w5Akey,w6Akey,w7Akey,w8Akey,w9Akey, labels = c("Fetal1", "Fetal2","Fetal3","Birth/Infan","Infan/Child","Child","Adolescent","Adult"),
          common.legend = TRUE, legend = "right")

ggsave(file="Sestan_PCA_WindowsAkey.pdf", width = 11.69, height = 8.27)
