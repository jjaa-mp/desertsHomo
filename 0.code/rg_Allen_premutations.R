```{r}
#Packages
library(biomaRt)
library(tidyr)
library(dplyr)
library("readxl")
library(ABAData)
library(ABAEnrichment)
library(regioneR)
library(GGally)
library(viridis)
library(ggpubr)
library(grid)
library(gridExtra)
library(lattice)
```
#Extracting gene names from coordinate regions via bioMart
```{r}
#Genes within Desert region coordinates
##Using hg19 genome
ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = 'hsapiens_gene_ensembl',
                        host = 'https://grch37.ensembl.org')
##Akey Deserts coordinates (Ensembl 1-based):
filterlist <- c("1:105400000:120600000", "3:74100000:89300000", "7:106200000:123200000","8:49400000:66500000")
##Select only protein-coding from Akey
resultsEnsembl=getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=filterlist,biotype="protein_coding"), mart = ensembl)
```

#Extracting gene expression data from Allen Brain Atlas
```{r}
##Loading dataset
data("dataset_adult") #Raül: dataset_adult (microarray data)
unique(dataset_adult$structure) #Checking structures

#(Perform enrichment on Genes from Deserts - Skip)
#input_hyper = data.frame(results$hgnc_symbol, is_candidate=1)
#res_devel = aba_enrich(input_hyper, dataset='5_stages')


#Selecting genes ID and structures present in dataset
id <- unique(dataset_adult$ensembl_gene_id)
st <- unique(dataset_adult$structure)
st_allen <- paste("Allen",st, sep=":") 
#Expression data - Getting raw data for all structures and genes
ab1 <- get_expression(structure_ids=st_allen, gene_ids = id, dataset='adult')
dataset_adult$
#Converting data to a dataframe with useful format for later
list1 = vector(mode="list")
abtemp<-t(ab1)


```
# NOW I'll perform the permutation test to evaluate the results
#First Performing the same process but now for 50 alternative sets of random regions

```{r}
#This mask will exclude the generation using samples from introgression deserts
mask<-data.frame(c("chr1","chr3","chr7","chr8"), c(105400000, 74100000,106200000,49400000), c(120600000, 89300000,123200000,66500000))
randReg50set1<-createRandomRegions(nregions = 50,length.mean = 15000000,length.sd = 1000000,genome = "hg19", mask = mask)
randReg50set2<-createRandomRegions(nregions = 50,length.mean = 15000000,length.sd = 1000000,genome = "hg19", mask = mask)
randReg50set3<-createRandomRegions(nregions = 50,length.mean = 15000000,length.sd = 1000000,genome = "hg19", mask = mask)
randReg50set4<-createRandomRegions(nregions = 50,length.mean = 15000000,length.sd = 1000000,genome = "hg19", mask = mask)

#Here I generate a proper dataframe that will be usefull later on
abtemp2 <- cbind(idGene = rownames(abtemp), abtemp)
rownames(abtemp2) <- 1:nrow(abtemp2)
abDeserts<-abtempFin<-abtemp2[abtemp2$idGene %in% resultsEnsembl$ensembl_gene_id,]
abDeserts$idGene<-NULL
drop(abDeserts$idGene)

# Initializing calculation of means per substructures with the proper dataframe
meansSubstruct<-colMeans(abDeserts)
rownames(meansSubstruct)
meansSubstructTemp <- cbind(idBrainStruct = rownames(meansSubstruct), meansSubstruct)
meansSubstructTemp <- cbind(idBrainStruct = rownames(meansSubstructTemp), meansSubstructTemp)
rownames(meansSubstructTemp) <- 1:nrow(meansSubstructTemp)


#This variable will contain all the ABAData per all the alternative permutations
randomABAData<-list()
for (i in 1:50){
  sprintf("%d randomset",i)
  desert1<-c(runValue(seqnames(randReg50set1)[i]),start(randReg50set1)[i],end(randReg50set1)[i])
  desert2<-c(runValue(seqnames(randReg50set2)[i]),start(randReg50set2)[i],end(randReg50set2)[i])
  desert3<-c(runValue(seqnames(randReg50set3)[i]),start(randReg50set3)[i],end(randReg50set3)[i])
  desert4<-c(runValue(seqnames(randReg50set4)[i]),start(randReg50set4)[i],end(randReg50set4)[i])
  filterlistTemp<-c(paste(desert1,collapse=":"),paste(desert2,collapse=":"),paste(desert3,collapse=":"),paste(desert4,collapse=":"))
  print(filterlistTemp)
  abtemp2 <- cbind(idGene = rownames(abtemp), abtemp)
  rownames(abtemp2) <- 1:nrow(abtemp2)
  resultsTemp=getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position","gene_biotype"),
                    filters = c("chromosomal_region","biotype"),
                    values = list(chromosomal_region=filterlistTemp,biotype="protein_coding"), mart = ensembl)
  abtempFin<-abtemp2[abtemp2$idGene %in% resultsTemp$ensembl_gene_id,]
  abtempFin$idGene<-NULL
  drop(abtempFin$idGene)
  #Here we safe all the data for each permutation
  randomABAData[[i]]<-abtempFin
  #Here we will safe the values of the means per substructures in a big dataframe
  meansSubStructRand<-colMeans(abtempFin)
  meansSubStructRandTemp <- cbind(idBrainStruct = rownames(meansSubStructRand), meansSubStructRand)
  meansSubStructRandTemp <- cbind(idBrainStruct = rownames(meansSubStructRandTemp), meansSubStructRandTemp)
  rownames(meansSubStructRandTemp) <- 1:nrow(meansSubStructRandTemp)
  meansSubstructTemp<-merge(meansSubstructTemp, meansSubStructRandTemp, by="idBrainStruct")
}
```

#Final calculations with the data from permutations

```{r}
#Here we calculate th t.test comparing the abDesert with each of the random samples
#In this case 50 t.test per each substructure.
brainstr<-DataFrame()
for (brainStruct in colnames(abDeserts)){
  suma<-0
  for (i in 1:50){
    val<-t.test(abDeserts[brainStruct],randomABAData[[i]][brainStruct],var.equal = TRUE)
    suma<-suma+val$p.value
    #tTests<-c(tTests,c(brainStruct),val)
  }
  brainstr[brainStruct]<-c(suma/50)
}

```
#Extracting the data to a proper csv
```{r}
lista <- as.list(as.data.frame(t(brainstr)))
brainStMean<-as.data.frame(t(brainstr))
brainStMean$group<-NULL
drop(brainStMean$group)

brainStMean<-brainStMean %>% 
  separate(group_name, c("Title", "id"),":")
brainStMean$Title<-NULL
drop(brainStMean$Title)

gene_exp_chr1<-DataFrame(read_excel("Data/genes_clean_desertic_regions_and_brain_regions.xlsx",sheet = "chr1"))
gene_exp_chr1<-gene_exp_chr1[c("id", "name", "topStructure")]
brainStMean<-merge(brainStMean,gene_exp_chr1,by=)
write.csv(brainStMean, "brainStMeanPerm50.csv", row.names = FALSE)
```

