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


- **Enrichment in Cerebrum and Cerebellum** We investigated the expression patterns across organs. We observe that genes present in deserts and expressed in the Cerebrum, Cerebellum, and Eye (neuroectoderm-derived structures) are significantly different from the rest of structures, while no structures passed the significance threshold for genes present in AkeyPey.
  + Mean expression of genes in Akey higher in Cerebrum, Cerebellum and Eye; bigger difference when considering genes in Akey and Pey in Cerebrum and Cerebellum (stats pending).
  + *Cell type enrichment*. No significant results were found for any cell type in neither organ (Cerebrum, Cerebellum, Eye).
- **414 structures** (in process).
- **ABA temporal progression**  Previously, Mafessoni2020 have shown that genes carrying Neanderthal-derived changes and expressed in the striatum during adolescence present a higher McDonald-Kreitman ratio. These same genes are over-represented in genomic regions were modern humans are significantly devoid of Neanderthal introgressed haplotypes (Maffesoni2020) [do they give the region coordinates?]. In addition, Vernot2016 [using a different range of introgressed regions] emphasized that genes within deserts are significantly enriched in the developing cerebral cortex and in the striatum at adolescence and adult stages. We investigated the temporal progression of expression patterns of genes within deserts of introgression and putative positively-selected regions in different brain structures accessing RNA-seq data from the Allen Brain Atlas (via the ABAData package; see Methods). We found that genes within deserts have a higher mean expression profile in the striatum during adolescence and, prominently, during infancy (stats pending; Figure). When considering genes present in both deserts and putative positively-selected regions, we found a clearly distinct profile for the cerebellum after prenatal stages, from infancy to adulthood (stats pending; Figure), while the rest of structures approximate a L-shaped trajectory. [also compare with raw data profile, clearly distinct]


## Stats

## Anova
Results on a series of aligned rank transform anovas (non-parametric) on the ABA temporal progression data can be seen [here](https://github.com/jjaa-mp/raul_tesina/tree/master/0.code/R_scripts/output/anova). It includes pairwise comparison in:
- Brain regions only
- Developmental stages only
- The interaction between brain region and developmental stages
For both the deserts of introgression and deserts + Peyrégne et al.

## Kruskal Wallis
A series of Kruskal Wallis tests show that the difference in expression values between brain regions (regardless of stages) in both files is not significant (p-values 0.7746, 0.9965). A post-hoc Dunn test also fails to highlight any region as well. A Wilcoxon signed-rank focusing on those stages significant in the anovas confirms this when ran specifically on the cerebellum and the prenatal-child stages explicitely.
