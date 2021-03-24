#Supplementary Tables at https://github.com/jjaa-mp/desertsHomo/tree/master/1.data/stat_freq
##Frequency data retrieved with permissions from Kuhlwilm M, Boeckx C. Sci Rep. 2019. http://www.nature.com/articles/s41598-019-44877-x

SummarySection1 <-read.csv("~/desertsHomo/1.data/stat_freq/Supplementary_Table_1.csv",sep = "\t")
dfreq <- read.csv("~/desertsHomo/1.data/stat_freq/Supplementary_Table_2.csv")

pairwise.wilcox.test(dfreq$freq, dfreq$group, p.adjust.method = "BH")
a <- pairwise.wilcox.test(dfreq$freq, dfreq$group, p.adjust.method = "BH")
a$p.value

##Total numbers for Chi-sq test
d <- as.table(rbind(c(591,3652-591), c(2587,31097-2587),  c(110,238-110)))
dimnames(d) <- list(group = c("Deserts","NoDeserts_chr", "Deserts+PosSel"), fixation = c("Fixed","No Fixed"))

chisq.test(d)

##Pairwise comparisons
p1 <- as.table(rbind(c(591,3652-591),  c(2587,31097-2587)))
dimnames(p1) <- list(group = c("Deserts","NoDeserts_chr"), fixation = c("Fixed","No Fixed"))
chisq.test(p1)

p2 <- as.table(rbind(c(591,3652-591), c(110,238-110)))
dimnames(p2) <- list(group = c("Deserts","Deserts+PosSel"), fixation = c("Fixed","No Fixed"))
chisq.test(p2)

p3 <- as.table(rbind( c(2587,31097-2587), c(110,238-110)))
dimnames(p3) <- list(group = c("NoDeserts_chr","Deserts+PosSel"), fixation = c("Fixed","No Fixed"))
chisq.test(p3)

pvals <- c(chisq.test(p1)$p.value,chisq.test(p2)$p.value,chisq.test(p3)$p.value)

p.adjust(pvals, method ="BH", n = length(pvals))


#Summary table as Figure
colnames(SummarySection1) <- c("Region", "Mean frequency", "% Fixed alleles")
SummarySection1 <- SummarySection1[c(2,1,3),]
levels(SummarySection1$Region) <- c(levels(SummarySection1$Region), "Non desertic")
SummarySection1$Region[SummarySection1$Region == 'No deserts'] <- 'Non desertic'
SummarySection1[,-1] <- round(SummarySection1[,-1], 4)

tabl_1 <- tableGrob(SummarySection1, rows=NULL)

#grid.arrange(tabl_1,tabl_1, ncol=2) #Share with circosplot



#HF strict filtering
SummarySection_strict <-read.csv("~/desertsHomo/1.data/stat_freq/ST_Summary_HFstrict_summary.csv",sep = "\t")

dst <- read.csv("~/desertsHomo/1.data/stat_freq/ST_HFstrict_frequency.csv", sep = ",")
pairwise.wilcox.test(dst$freq, dst$group, p.adjust.method = "BH")
a1<-pairwise.wilcox.test(dst$freq, dst$group, p.adjust.method = "BH")
a1$p.value
##Total numbers for Chi-sq test
d <- as.table(rbind(c(587,2246-587), c(2544,17783-2544),  c(110,229-110)))
dimnames(d) <- list(group = c("Deserts","NoDeserts_chr", "Deserts+PosSel"), fixation = c("Fixed","No Fixed"))

chisq.test(d)

##Pairwise comparisons
p1 <- as.table(rbind(c(587,2246-587),   c(2544,17783-2544)))
dimnames(p1) <- list(group = c("Deserts","NoDeserts_chr"), fixation = c("Fixed","No Fixed"))
chisq.test(p1)

p2 <- as.table(rbind(c(587,2246-587), c(110,229-110)))
dimnames(p2) <- list(group = c("Deserts","Deserts+PosSel"), fixation = c("Fixed","No Fixed"))
chisq.test(p2)

p3 <- as.table(rbind(  c(2544,17783-2544),c(110,229-110)))
dimnames(p3) <- list(group = c("NoDeserts_chr","Deserts+PosSel"), fixation = c("Fixed","No Fixed"))
chisq.test(p3)

pvals <- c(chisq.test(p1)$p.value,chisq.test(p2)$p.value,chisq.test(p3)$p.value)

p.adjust(pvals, method ="BH", n = length(pvals))