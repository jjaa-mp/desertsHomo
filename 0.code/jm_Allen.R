```{r}
#Packages
library(biomaRt)
library(ABAData)
library(ABAEnrichment)
library(dplyr)

#Genes within Desert region coordinates
##Using hg19 genome
ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                        dataset = 'hsapiens_gene_ensembl',
                        host = 'https://grch37.ensembl.org')
##Akey Deserts coordinates (Ensembl 1-based):
filterlist <- c("1:105400000:120600000", "3:74100000:89300000", "7:106200000:123200000","8:49400000:66500000")
##Select only protein-coding from Akey
results=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=filterlist,biotype="protein_coding"), mart = ensembl)
##Save
write.csv(results, file="2020_Genes_in_Deserts.csv")

##Pey coordinates (Ensembl 1-based):  ##CHECK COORDINATES DEPENDING ON FILE (0-BASED vs 1-BASED)
pey_coords <- read.delim("~/2020_pey_coords.bed", header=FALSE)
###Preparing bed file for input in bioMart
pey_coords$V1 <- gsub("chr", "\\1", pey_coords$V1)
pey_coords
df <- paste(pey_coords$V1, pey_coords$V2, pey_coords$V3, sep = ":")
results_pey=getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position","gene_biotype"),
              filters = c("chromosomal_region","biotype"),
              values = list(chromosomal_region=df,biotype="protein_coding"), mart = ensembl)
###Genes within Akey and Pey - RESULT
both <- results_pey[results_pey$hgnc_symbol %in% results$hgnc_symbol,]
###Cleaning - To be used later
both <-  both[!(is.na(both$hgnc_symbol) | both$hgnc_symbol==""), ]


#Extracting gene expression data from Allen Brain Atlas

##Loading dataset
data("dataset_5_stages") #Raül: dataset_adult (microarray data)
unique(dataset_5_stages$structure) #Checking structures

#(Perform enrichment on Genes from Deserts - Skip)
#input_hyper = data.frame(results$hgnc_symbol, is_candidate=1)
#res_devel = aba_enrich(input_hyper, dataset='5_stages')


#Selecting genes ID and structures present in dataset
id <- unique(dataset_5_stages$ensembl_gene_id)
st <- unique(dataset_5_stages$structure)
st_allen <- paste("Allen",st, sep=":") 
#Expression data - Getting raw data for all structures and genes
ab <- get_expression(structure_ids=st_allen, gene_ids = id, dataset='5_stages')

#Converting data to a dataframe with useful format for later
list1 = vector(mode="list")
for (r in 1:length(ab)){
  ab[[r]] <- t(ab[[r]]) #transpose
  list1 <- get_name(colnames(ab[[r]])) #change Allen:XXXX to e.g. M1C_primary motor cortex, etc
  colnames(ab[[r]]) <- list1
  ab[[r]] <- as.data.frame(ab[[r]])
}


#Selection of >q75 to select genes with high expression
q75 = vector(mode="list", length = length(ab))
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    q75[[i]][[h]] <- ab[[i]][h] %>% filter(ab[[i]][h] > quantile(ab[[i]][[h]], 0.75))
  }
}
#Changing row names ENSG to hgnc_symbol (gene names; via bioMart); ordering values based on expression
##ab is a list of 5 elements (i); each element or 'sublist' has 16 dataframes (h)
##q75[[1]][[3]][[1]] would denote  primer sublista (1), tercer dataframe (3), columna 1 (0-based)
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(q75[[i]][[h]]),mart=ensembl)
    q75[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(q75[[i]][[h]]), G_list$ensembl_gene_id)]
    q75[[i]][[h]] <-q75[[i]][[h]][order(q75[[i]][[h]][[1]], decreasing = TRUE), ]
    
  }
}


#q9 - Skip
#q9 = vector(mode="list", length = length(ab))
#for (i in 1:length(ab)){
  #for (h in 1:length(names(ab[[i]]))){
    #q9[[i]][[h]] <- ab[[i]][h] %>% filter(ab[[i]][h] > quantile(ab[[i]][[h]], 0.9))
  #}
#}

#library('biomaRt')
#mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#for (i in 1:length(ab)){
  #for (h in 1:length(names(ab[[i]]))){
    #G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(q9[[i]][[h]]),mart=ensembl)
    #q9[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(q9[[i]][[h]]), G_list$ensembl_gene_id)]
    #q9[[i]][[h]] <-q9[[i]][[h]][order(q9[[i]][[h]][[1]], decreasing = TRUE), ]
  #}
#}

#qlower
#qlower = vector(mode="list", length = length(ab))
#for (i in 1:length(ab)){
  #for (h in 1:length(names(ab[[i]]))){
    #qlower[[i]][[h]] <- ab[[i]][h] %>% filter(ab[[i]][h] < quantile(ab[[i]][[h]], 0.25))
  #}
#}
#for (i in 1:length(ab)){
  #for (h in 1:length(names(ab[[i]]))){
    #G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"),values=rownames(qlower[[i]][[h]]),mart=ensembl)
    #qlower[[i]][[h]]['gene_name'] <-  G_list$hgnc_symbol[match(rownames(qlower[[i]][[h]]), G_list$ensembl_gene_id)]
    #qlower[[i]][[h]] <-qlower[[i]][[h]][order(qlower[[i]][[h]][[1]], decreasing = TRUE), ]
  #}
#}

#Intersecting q75 (high expression) with Akey
akey075 = vector(mode="list", length = length(ab)) #Creating empty list with 5 elements as ab
#(h in 1:length(names(ab[[i]]))) to generate same number of dataframes (in this case 16) as original in ab 
for (i in 1:length(ab)){
  for (h in 1:length(names(ab[[i]]))){
    akey075[[i]][[h]] <- q75[[i]][[h]][q75[[i]][[h]][[2]] %in% results$hgnc_symbol,] #in Akey
    akey075[[i]][[h]] <-  akey075[[i]][[h]][!(is.na(akey075[[i]][[h]][[2]]) | akey075[[i]][[h]][[2]]==""), ] #Cleaning
  }
}

##Akey09
#akey = vector(mode="list", length = length(ab))
#for (i in 1:length(ab)){
  #for (h in 1:length(names(ab[[i]]))){
    #akey[[i]][[h]] <- q9[[i]][[h]][q9[[i]][[h]][[2]] %in% results$hgnc_symbol,]
    #akey[[i]][[h]] <-  akey[[i]][[h]][!(is.na( akey[[i]][[h]][[2]]) | akey[[i]][[h]][[2]]==""), #]
  #}
#}
##AkeyLower
#AkeyLower = vector(mode="list", length = length(ab))
#for (i in 1:length(ab)){
  #for (h in 1:length(names(ab[[i]]))){
    #AkeyLower[[i]][[h]] <- qlower[[i]][[h]][qlower[[i]][[h]][[2]] %in% results$hgnc_symbol,]
    #AkeyLower[[i]][[h]] <-  AkeyLower[[i]][[h]][!(is.na( AkeyLower[[i]][[h]][[2]]) |  AkeyLower[[i]][[h]][[2]]==""), ]
  #}
#}

#Save results - Genes in deserts and pos. selection with high expression
for (i in 1:length(ab)){
  for (h in 1:length(bothq75[[i]])){
    if (dim(bothq75[[i]][[h]])[1] == 0) next #To skip empty dataframes
    write.xlsx(bothq75[[i]][[h]], file="ABA_preliminar_AkeyPey_highExpr.xlsx", sheetName=paste(toString(i),toString(names(bothq75[[i]][[h]][1])), sep="_"), append = TRUE) #Add age data at the beginning of the sheet name
  }
}


#--SKIP--
#Plot example
##Highlight specific gene
##Prenatal - WDR47
library(dplyr)
library(readxl)
library(ggplot2)
DLPFC_pre <- read_excel("~/ABA_preliminar.xlsx", sheet = 2)
highlight1 <-  as.data.frame(DLPFC_pre %>% filter(DLPFC_pre$gene_name == "WDR47"))
ggplot(DLPFC_pre, aes(x=DLPFC_pre$gene_name, y = DLPFC_pre$`DFC_dorsolateral prefrontal cortex`)) + 
  geom_point(binaxis='y', stackdir='center',stackratio=1.5)+theme(axis.text.x = element_text(angle=90), legend.position = "none") +theme(axis.text.x = element_text(margin = margin(t = 4, r = 20, b = 0, l = 0))) +
  geom_point(data=highlight1, aes(x=highlight1$gene_name, y = highlight1$`DFC_dorsolateral prefrontal cortex`, colour = "red", size=3))+
  labs(x="Gene name",y="expression", title ="Allen Brain Atlas - Prenatal DLPFC - WDR47")
# + coord_flip()

```

#Table Cell Atlas - Science
```{sh}
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156793/suppl/GSE156793%5FS8%5FDE%5Fgene%5Fcells%2Ecsv%2Egz 
```

```{r}
library(dplyr)
df <- read.csv(file = "~/GSE156793_S8_DE_gene_cells.csv.gz")
df1 <- df[which(df$organ=='Cerebellum'), ]
write.csv(df1, file="CellAtlas_GSE156793_cerebellum.csv", row.names = FALSE)

df1$gene_short_name <- gsub("\\'", "", df1$gene_short_name)
df1subset <- df1[df1$gene_short_name %in% results$hgnc_symbol,]
df1subset <- df1subset[order(df1subset$max.cluster),]

df1subsetboth <- df1[df1$gene_short_name %in% both$hgnc_symbol,]
df1subsetboth <- df1subsetboth[order(df1subsetboth$max.cluster),]
df1subset <- df1subset[order(df1subset$max.cluster),]


write.csv(df1subset, file="CellAtlas_GSE156793_inAkey.csv", row.names = FALSE)
write.csv(df1subsetboth, file="CellAtlas_GSE156793_inAkeyPey.csv", row.names = FALSE)


```


