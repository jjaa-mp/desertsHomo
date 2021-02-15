#Supplementary Table 1
##Frequency data retrieved with permissions from Kuhlwilm M, Boeckx C. Sci Rep. 2019. http://www.nature.com/articles/s41598-019-44877-x

SummarySection1 <-read.csv("Supplementary_Table_1.csv",sep = "\t")
dfreq <- read.csv("Supplementary_Table_2.csv")

pairwise.wilcox.test(dfreq$freq, dfreq$group, p.adjust.method = "BH")


##Total numbers for Chi-sq test
d <- as.table(rbind(c(591,3652-591), c(2587,31097-2587),  c(110,238-110)))
dimnames(d) <- list(group = c("Deserts","NoDeserts_chr", "Deserts+PosSel"), fixation = c("Fixed","No Fixed"))
d

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
