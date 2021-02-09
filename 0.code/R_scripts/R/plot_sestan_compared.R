plot_sestan_compared <- function(sestan_a, sestan_apey){
  raw <- read.csv("file_dependencies/df_raw.csv")
  raw <- raw[2:10] #see cluster script
  
 
  colnames(raw) <- colnames(sestan_a) 
  #silly change in cluster of colnames
  raw$dataset <- c("raw")
  sestan_a$dataset <- c("akey")
  sestan_apey$dataset <- c("akeypey")
  
  tot_pl <- base::rbind(sestan_a, sestan_apey, raw)
  tot_pl <- as_tibble(tot_pl)
  tot_pl <- tot_pl %>% 
    dplyr::mutate(dataset=as.character(dataset))
  
  levels(tot_pl$dataset) <-  c("akey", "akeypey", "raw")
  
  n <- ggparcoord(tot_pl,
                 columns = 2:8, 
                 groupColumn = 1, 
                 showPoints = TRUE, 
                 scale = "globalminmax",
                 title="Structures ABA") +
    theme(plot.title = element_text(size=10), legend.position = "none", axis.text.x = element_text(angle = 45,hjust = 1)) +
    xlab("") +
    ylab("expression")
  
  
  
  n <- ggparcoord(tot_pl,
                 columns = 2:8, 
                 groupColumn = 1,
                 showPoints = TRUE, 
                 scale = "globalminmax",
                 title="Structures Sestan", 
                 mapping=aes(color=as.factor(dataset))) +
    theme(plot.title = element_text(size=10), axis.text.x = element_text(angle = 45,hjust = 1)) +
    xlab("") +
    ylab("expression") +
    labs(color="Dataset")
  
  n + facet_wrap(~Structure)+
    scale_color_discrete(name="Dataset",
                         labels=unique(tot_pl$dataset))
  
  an <- n + facet_wrap(~Structure)+
    scale_color_discrete(name="Dataset",
                        labels=unique(tot_pl$dataset))
  ggsave(file="output/Sestan_temporal_Structures.pdf", an, width = 11.69, height = 8.27, units = "in")
}