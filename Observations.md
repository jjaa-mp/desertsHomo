# Observations
## Highlights from the data analysis
Key for me: make clear why we focus on Akey, why we further focus on Pey (I suggest we just focus on the 'extended', not 'core' set), and point to differences from the early, tentative Vernot conclusion

For intro: go back to work by Martin (his FOXP2 paper, and his Bioessays piece) to introduce the relevance of introgression deserts fo human uniqueness (in this content, see also Evo Letters 2020 piece on "genomic islands of divergence/speciation" in introgression contexts). Also briefly recap Vernot results, then say why we focus on Akey, and further examine the overlap with Pey. Striatum/CBL, like Gunz/Tilot?

From results:

- **SNP count and frequency** (*State importance of MartinBoeckx approach*). We evaluated the mean frequency of the HF variants identified in deserts and putative positively-selected regions. Being the differences small, we found significant results when comparing the mean frequency of variants within deserts and non-desertic regions, as well as to those falling within both deserts and positive selected regions -where we find a very high degree of fixation (>0.99)- (Wilcoxon rank sum test with Benjamini-Hochberg multiple comparison correction). Additionally, we found that the ratio of fixed alleles to total alleles is much higher in both deserts and positive selected regions (~46%) than in deserts (~16%) and non-desertic (~8%) regions (being all these differences statistically significant; Chi-squared test, p-value adjusted by Benjamini-Hochberg correction).

| Regions  | Frequency | % Fixed alleles |
| ------------- | ------------- | ------------- |
| Deserts  | 0.969864  |  16.182913  |
| No deserts (chr)  | 0.963944  | 8.319130  |
| Deserts + Pey  | 0.994599  | 46.218487  |


- *Cell Atlas* 
	+ *Cell type enrichment*. No significant results were found for any cell type in neither organ (Cerebrum, Cerebellum, Eye).

- **414 structures** [Need to be udpated] PCA performed on the whole dataset revealed cerebellum as a clearly distinct cluster (inside it, right/left separation was noticed too), as well as broad brain cortical and nuclei clusters. Permutations to compare expression levels between Chen, Chen+Pey in progress.


- **Sestan permutations**

- **Sestan temporal progression**  [Previously, Mafessoni2020 have shown that genes carrying Neanderthal-derived changes and expressed in the striatum during adolescence present a higher McDonald-Kreitman ratio. These same genes are over-represented in genomic regions were modern humans are significantly devoid of Neanderthal introgressed haplotypes (Maffesoni2020) [do they give the region coordinates?]. In addition, Vernot2016 [using a different range of introgressed regions] emphasized that genes within deserts are significantly enriched in the developing cerebral cortex and in the striatum at adolescence and adult stages. We investigated the temporal progression of expression patterns of genes within deserts of introgression and putative positively-selected regions in different brain structures accessing RNA-seq data from Sestan.] 
  + **Striatum**. We found that genes within deserts present highest median expression during childhood, largely decreasing afterwards. We notice, however, that genes within deserts present higher expression profile for the cerebellum prenatally, and for the MD thalamus during and after childhood. When subsetting for genes within both deserts and positive selection, the striatum mostly recapitulates the cup-shaped trajectory found in the global profile of this structure.
  + **Cerebellum**. We found that genes within deserts alone, and opposed to the global expression profile, present a median expression higher at prenatal stages (fetal stages 1,2,3). When considering genes present in both deserts and putative positively-selected regions, the expression profile increases as early as fetal stage 3 until it reaches the highest median expression of all structures from childhood to adulthood. 
  + **Akey & AkeyPey - General comment**. When comparing the expression profiles of genes in these two different subsets we observe that, in general, the trajectories of genes in both AkeyPey are more variable peaking at prenatal stages and decreasing later on, with the exception of the cerebellum. 
  + **PC & Pairwise Euclidean distances**.


---
## Summary of Methods

*Degree of fixation*. To evaluate the degree of fixation of Homo sapiens-specific variants in regions of interest, we processed a dataset containing SNCs derived in the modern human lineage at fixed or nearly fixed frequency (above 90% in present-day human populations) where high-quality Neanderthal/Denisovan genomes present the ancestral allele (from MartinBoeckx), and integrated it with genomic data of regions depleted of archaic haplotypes (from Chen2020) and putative positively-selected regions (from Pey2017). Statistical evaluation of differences in mean frequency and percentage of fixed alleles  was performed with a Wilcoxon rank sum test with and Chi-squared test, respectively, with Benjamini-Hochberg correction for multiple comparisons.


*Sestan. Temporal progression.* Genes within deserts of introgression or putative positively-selected regions were obtained via BioMart R package, inputing the respective genomic region coordinates (extracted from Chen2020 and Pey2017) and filtering by protein-coding genes. Publicly available transcriptomic data of the human brain at different developmental stages was retrieved from (Sestan, psychencode). RPKM normalized counts were log-transformed and then subsetted to select genes either in deserts of introgression or in both deserts and putative positively-selected regions. The complete log-transformed, RPKM normalized count matrix was subsetted to select genes with median expression value higher than 2, following (Sestan).

---


---
---
**Stats run on ABA data** 

**Non-parametric stats on ABA temporal progression**

These stats were run with the raw data, and don't require normal data (since they are non-parametric). However, using log values and adding an expression cutoff might change them if we decide to do so.
	
  + Anova: Results on a series of aligned rank transform anovas (non-parametric) on the ABA temporal progression data can be seen [here](https://github.com/jjaa-mp/raul_tesina/tree/master/0.code/R_scripts/output/anova). It includes pairwise comparison in:
		- Brain regions only
		- Developmental stages only
		- The interaction between brain region and developmental stages
	(For both the deserts of introgression and deserts + Peyrégne et al.)

- **Does cerebellum have a particular profile in ABA data?**

Motivated by the peculiar shape of the cerebellum in the means plots, we did a series of a series of Kruskal Wallis tests show that the difference in expression values between brain regions (regardless of stages) in both files is not significant (p-values 0.7746, 0.9965). A post-hoc Dunn test also fails to highlight any region as well. A Wilcoxon signed-rank focusing on those stages significant in the anovas confirms this when ran specifically on the cerebellum and the prenatal-child stages explicitely. However, account for genes as groups does solve this, as a Friedmann test made with the data of each stage independently of Chen and Chen + Pey shows significance for everything except the prenatal stage of Chen + Pey (corrected by Bonferroni). 

