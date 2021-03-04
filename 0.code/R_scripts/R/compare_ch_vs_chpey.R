compare_ch_vs_chpey <- function(sestan){
  ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'https://grch37.ensembl.org')
  
  Deserts <- c("1:105400000-120600000", "3:74100000-89300000",
               "7:106200000-123200000","8:49400000-66500000")
  
  DesertsPosSel <- c("1:113427676-113560554", "1:114641362-114645248",
                     "1:119322276-119387279", "3:77027847-77034264",
                     "7:106877730-107233808", "7:116762909-116773234",
                     "7:120147456-120174406", "7:122320035-122406480") #can be obtained from .bed file in circos plot files, if you wish so
  
  
  genesdeserts <- getBM(attributes = c("hgnc_symbol"),
                        filters = c("chromosomal_region","biotype"),
                        values = list(Deserts,biotype="protein_coding"), mart = ensembl)
  
  genesdesertspos <- getBM(attributes = c("hgnc_symbol"),
                           filters = c("chromosomal_region","biotype"),
                           values = list(DesertsPosSel,biotype="protein_coding"), mart = ensembl)
  
  
  #Now, determining desert specific genes NOT under possel:
  genesnotin_pos <- anti_join(genesdeserts, genesdesertspos)
  
  #Select random 488 genes in positive selection + deserts
  # 488 because that's the body of genes in deserts not under possel
  set.seed(1)
  genesRandPS <- sample(genesdesertspos$hgnc_symbol, 488)
  genesRandPS_expr <- sestan %>% 
    filter(gene_name %in% genesRandPS) #466 instead of 488 due to name divergence
  
  #Cutoff + normalization:
  genesRandPS_expr <- melt(genesRandPS_expr)
  genesRandPS_expr <- genesRandPS_expr %>% 
    filter(value > 2) %>% 
    dplyr::mutate(value = log2(value))
  
  genesRandPS_expr <- separate(genesRandPS_expr, variable, c("stage", "structure"), sep = "[.]")
  genesRandPS_expr$stage <- as.factor(genesRandPS_expr$stage)
  
  # Now for those NOT under positive selection:
  genes_notunder_ps <- sestan %>% 
    filter(gene_name %in% genesnotin_pos$hgnc_symbol) #469 instead of 488
  
  genes_notunder_ps <- melt(genes_notunder_ps)
  #Cutoff + normalization (but for the other dataset)
  genes_notunder_ps <- genes_notunder_ps %>% 
    filter(value > 2) %>% 
    dplyr::mutate(value = log2(value))
  genes_notunder_ps <- separate(genes_notunder_ps, variable, c("stage", "structure"), sep = "[.]")
  
  p <- NULL
  p$'Deserts (excluding positive selection)'  <- genes_notunder_ps$value
  p$'Deserts and positive selection' <- genesRandPS_expr$value
  p <- melt(p)
  box <- ggplot(p, aes(L1, value)) +
    theme_minimal() +
    labs(y = "Log2 expression per gene (n=466)", x= "Category") +
    geom_boxplot()
  ggsave("boxplot.pdf", box, width = 8, height = 8)
  #stats
  kruskal.test(value~L1, data = p)
  
  
  p <- NULL
  p <- genes_notunder_ps[3:4]
  p <- rbind(p, genesRandPS_expr[3:4])
  datasource1 <- rep("Not under ps", length(genes_notunder_ps$value))
  datasource2 <- rep("Under ps", length(genesRandPS_expr$value))
  p$datasource <- c(datasource1, datasource2)
  
  violins <- ggplot(p, aes(datasource, value, fill = datasource)) +
    theme_minimal() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Log2 expression per gene (n=466)", x= "Category") +
    geom_violin() + 
    facet_wrap(vars(structure),  ncol = 6)
  ggsave("violin_plots.pdf", violins, width = 8, height = 8)
  
  #Stats!
  kruskal.test(value~structure, data = p) #first approximation: structure regardless of datasource significantly different
  
  model <-  lme(value ~ datasource, random=~1|structure,
                data=p,
                method="REML")
  
  anova.lme(model,
            type="sequential",
            adjustSigma = FALSE) # Data source is significant
  
  leastsquare <- lsmeans(model, pairwise ~ datasource,
                         adjust="tukey")
  leastsquare  # agreement
  
  # Test significance of random variable
  model.fixed <- gls(value ~ datasource,
                     data=p,
                     method="REML")
  
  anova(model, model.fixed) #0.042, meaning difference between structures is significant
  # (independently of where in the deserts they are)
  
  
}


#taking into account stages
p <- NULL
p <- genes_notunder_ps[2:4]
p <- rbind(p, genesRandPS_expr[2:4])
