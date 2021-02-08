library(tidyverse)
library(GGally)
library(gridExtra)
library(ggpubr)
library(grid)

df_raw <- read_csv("df_raw.csv")

colnames(df_raw) <- c("Index", "Structure", "Fetal1", "Fetal2", "Fetal3", "Birth/Infancy", "Infancy/Child", "Child", "Adolescence", "Adult")

levels(colnames(df_raw)) <- c("Index", "Structure", "Fetal1", "Fetal2", "Fetal3", "Birth/Infancy", "Infancy/Child", "Child", "Adolescence", "Adult")



a <- ggparcoord(df_raw,
               columns = 3:10, groupColumn = 2, showPoints = TRUE, 
               scale = "globalminmax",
               title="Genes in Deserts and Pos Sel- VFC (green) & AMY (black)")+
  theme_minimal() +
  scale_color_manual(values = c( "#ABABAB", "#000000", 
                                 "#ABABAB", "#ABABAB",
                                 "#ABABAB", "#ABABAB", 
                                 "#ABABAB",  "#ABABAB", 
                                 "#ABABAB", "#ABABAB",
                                 "#ABABAB", "#ABABAB", 
                                 "#ABABAB", "#ABABAB", 
                                 "#ABABAB", "#238b45")) +
  theme(plot.title = element_text(size=10),legend.position = "none", axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") +
  ylab("Expression")
b <- ggparcoord(df_raw,
                columns = 3:10, groupColumn = 2, showPoints = TRUE,
                scale = "globalminmax",title="Striatum")+
  theme_minimal() +
  scale_color_manual(values = c( "#ABABAB", "#ABABAB", 
                                 "#ABABAB", "#ABABAB",
                                 "#ABABAB", "#ABABAB", 
                                 "#ABABAB",  "#ABABAB", 
                                 "#ABABAB", "#ABABAB",
                                 "#ABABAB", "#ABABAB", 
                                 "#ABABAB", "#FF0000", 
                                 "#ABABAB", "#ABABAB"))+
  theme(plot.title = element_text(size=10),legend.position = "none", axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") +
  ylab("Expression")
c<-ggparcoord(df_raw,
              columns = 3:10, groupColumn = 2, showPoints = TRUE, scale = "globalminmax",
              title="Cerebellum") +
  theme_minimal() +
  scale_color_manual(values = c( "#ABABAB", "#ABABAB", 
                                 "#0000FF", "#ABABAB",
                                 "#ABABAB", "#ABABAB", 
                                 "#ABABAB",  "#ABABAB",
                                 "#ABABAB", "#ABABAB",
                                 "#ABABAB", "#ABABAB", 
                                 "#ABABAB", "#ABABAB", 
                                 "#ABABAB", "#ABABAB"))+
  theme(plot.title = element_text(size=10),legend.position = "none", axis.text.x = element_text(angle = 45,hjust = 1)) +
  xlab("") +
  ylab("Expression")

a1<-arrangeGrob(a, left=textGrob("A"))
b1<-arrangeGrob(b, left =textGrob("B"))
c1<-arrangeGrob(c, left=textGrob("C"))
grid.arrange(a, arrangeGrob(b, c), ncol = 2)

plot <-grid.arrange(a1, b1, c1, ncol = 2, layout_matrix = rbind(c(1, 1, 2), c(1, 1, 3)))
ggsave("raw_sestan.pdf", plot, width = 8, height = 8)
