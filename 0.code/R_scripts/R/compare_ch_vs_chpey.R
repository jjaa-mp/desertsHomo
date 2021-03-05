compare_ch_vs_chpey <- function(sestan){
  ensembl <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', 
                     dataset = 'hsapiens_gene_ensembl',
                     host = 'https://grch37.ensembl.org')
  
  Deserts <- c("1:105400000:120600000", "3:74100000:89300000",
               "7:106200000:123200000","8:49400000:66500000")
  
  DesertsPosSel <- c("1:113427676:113560554", "1:114641362:114645248",
                     "1:119322276:119387279", "3:77027847:77034264",
                     "7:106877730:107233808", "7:116762909:116773234",
                     "7:120147456:120174406", "7:122320035:122406480") #can be obtained from .bed file in circos plot files, if you wish so
  
  
  genesdeserts <- getBM(attributes = c("hgnc_symbol"),
                        filters = c("chromosomal_region","biotype"),
                        values = list(Deserts,biotype="protein_coding"), mart = ensembl)
  
  genesdesertspos <- getBM(attributes = c("hgnc_symbol"),
                           filters = c("chromosomal_region","biotype"),
                           values = list(DesertsPosSel,biotype="protein_coding"), mart = ensembl)
  
  
  #Now, determining desert specific genes NOT under possel:
  genesnotin_pos <- anti_join(genesdeserts, genesdesertspos)
  
  genesdes_expr <- sestan %>% 
    filter(gene_name %in% genesnotin_pos$hgnc_symbol) 
  
  #Cutoff + normalization:
  genesdes_expr <- melt(genesdes_expr)
  genesdes_expr <- genesdes_expr %>% 
    filter(value > 2) %>% 
    dplyr::mutate(value = log2(value))
  
  genesdes_expr <- separate(genesdes_expr, variable, c("stage", "structure"), sep = "[.]")
  genesdes_expr$stage <- as.factor(genesdes_expr$stage)
  
  # Now for those under positive selection:
  genesdesertspos <- sestan %>% 
    filter(gene_name %in% genesdesertspos$hgnc_symbol) #469 instead of 488
  
  genesdesertspos <- melt(genesdesertspos)
  #Cutoff + normalization (but for the other dataset)
  genesdesertspos <- genesdesertspos %>% 
    filter(value > 2) %>% 
    dplyr::mutate(value = log2(value))
  genesdesertspos <- separate(genesdesertspos, variable, c("stage", "structure"), sep = "[.]")
  
  p <- NULL
  p$'Genes in deserts (no positive selection)'  <- genesdes_expr$value
  p$'Deserts and positive selection' <- genesdesertspos$value
  p <- melt(p)
  violin1 <- ggplot(p, aes(value, L1)) + 
    theme_minimal() +
    geom_violin()
  ggsave("violin_basic.pdf", violin1, width = 8, height = 8)
  #stats
  kruskal.test(value~L1, data = p)
  
  
  p <- NULL
  p <- genesdes_expr[3:4]
  p <- rbind(p, genesdesertspos[3:4])
  datasource1 <- rep("Not under ps", length(genesdes_expr$value))
  datasource2 <- rep("Under ps", length(genesdesertspos$value))
  p$datasource <- c(datasource1, datasource2)
  
  violins <- ggplot(p, aes(datasource, value, fill = datasource)) +
    theme_minimal() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    labs(y = "Log2 expression per gene (n=466)", x= "Category") +
    geom_violin() + 
    facet_wrap(vars(structure),  ncol = 6)
  ggsave("violin_structures.pdf", violins, width = 8, height = 8)
  
  #Stats!
  kruskal.test(value~structure, data = p) #Structure regardless of datasource not significantly different
  
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

  
  #taking into account stages
  p <- NULL
  p <- genesdes_expr[2:4]
  p <- rbind(p, genesdesertspos[2:4])
  datasource1 <- rep("Not under ps", length(genesdes_expr$value))
  datasource2 <- rep("Under ps", length(genesdesertspos$value))
  p$datasource <- c(datasource1, datasource2)
  p$stage <- p$stage %>% 
    stringr::str_replace_all("HSB153|HSB150|HSB113|HSB103|HSB149|HSB114", "fetal1" ) %>% 
    stringr::str_replace_all("HSB178|HSB154|HSB96|HSB97", "fetal2" ) %>% 
    stringr::str_replace_all("HSB98|HSB107|HSB92|HSB159", "fetal3" ) %>% 
    stringr::str_replace_all("HSB155|HSB194|HSB121|HSB132|HSB139", "birth_inf" ) %>% 
    stringr::str_replace_all("HSB131|HSB171|HSB122|HSB143|HSB173", "inf_child" ) %>% 
    stringr::str_replace_all("HSB172|HSB118|HSB141|HSB174|HSB175", "child" ) %>% 
    stringr::str_replace_all("HSB124|HSB119|HSB105|HSB127", "adolescence" ) %>%   
    stringr::str_replace_all("HSB130|HSB136|HSB126|HSB145|HSB123|HSB135", "adult" ) 
  
  p <- p %>% # taking out prenatal 
    filter(stage != ("HSB148")) %>% 
    filter(stage != ("HSB112"))
  
  violin_stages <- ggplot(p, aes(value, stage)) + 
    theme_minimal() +
    geom_violin()
  ggsave("violin_stages.pdf", violin_stages, width = 8, height = 8)
  
  model <-  lme(value ~ datasource, random=~1|stage,
                data=p,
                method="REML")
  
  anova.lme(model,
            type="sequential",
            adjustSigma = FALSE) # stage is significant
  
  leastsquare <- lsmeans(model, pairwise ~ datasource,
                         adjust="tukey")
  leastsquare  # agreement
  
  # Test significance of random variable
  model.fixed <- gls(value ~ datasource,
                     data=p,
                     method="REML")
  
  anova(model, model.fixed) # stage does play a role
  
  #But which stages?
  sig <- aov(value ~ stage*datasource, data = p) %>% tukey_hsd()
  write.csv(sig, "posthoc_stages_chvschpey.csv")
}


