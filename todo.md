
I included tasks here as questions we aim to resolve. I didn't specifically included the answers, but I marked them for whenever we have adressed them overall. 
Feel free to add more, or to add/modify the tasks. Also, I'm not explicitely including the answers to the solved tasks, to keep this clean in terms of logic.

# Deserts of introgression vs the rest of the genome
[X] are results from previous studies, namely Vernot, similar to what we get?
[X] Are deserts of introgression different from other regions in the genome?

# Deserts of introgression vs regions under positive selection within
[ ] Are regions under positive selection different in their brain expression compared to any other region in the deserts?
[ ] Is this difference noticeable in different brain regions?
[ ] Is it stage-sensitive?

# Regions unders positive selection VS those specifically in deserts
Not adressed so far

# Differences between deserts of introgression
[ ] Is chromosome 7 different in its expression profiles given previous literature on its human uniqueness?


# Differences between genes
[X] Visual representation within the deserts
[ ] What are the differences in the gene trajectories across regions? (stats)
[ ] Are some genes specifically more expressed in a certain brain ragion?
[ ] Have certain genes in the deserts a higher than expected mean expression? 

# Differences across stages of development in different sets
[X] Versus the whole genome
[ ] Between sets (deserts vs pos selection windows within?)
[X] Between genes

	+ AA: To what resolution do we want to answer this question?
	+ I propose it's interesting in terms of 

--- 

# Previous to do:

- [X] Permutation tests for ABA 414 structures dataset (Alejandro) (done)

- [X] **Sestan data**.
    + [X] Akey+Pey (Raül)
    + [X] Akey alone (Juan)
    + [X] Entire dataset (ubics cluster, Alejandro)
    
 - [X] **Sestan data** *specific plots* (Juan)
    + [X] Structure-by-structure plot
    + [X] Gene-by-gene (ChenPey) plot
  
- [X] **Stats** for ABA temporal progression (data updated with log2(rpkm+1)) (Juan) 
    + AkeyPey: All stages significant except 1; Akey: all stages

- [X] **Stats** for Sestan temporal progression (data with log2(rpkm+1))  (Juan)
    +Friedman + Conover (AkeyPey: significant differences at three stages; Akey: all stages; mult comp cannot spot particular structures)
    + [X] Sestan permutation runs (AA)
    
- [X] PCA for each stage (Raül)
    - [X] Pairwise distances from PC clusters (Raül)
    - [X] PCA for Sestan raw (Juan)
    - [X] Pairwise distances from PC clusters - raw data
    - [X] Stats on pairwise distances

  
- [ ] **Permutation subtasks** (AA)
    - [X] Difference plots (Sestan)
    - [X] Difference plots (ABA)
    - [X] Stats for diff plots
    - [X] Are some structures more prone to be at the top or bottom? Some kind of enrichment test w/ permutations / solved with outliers test (different kind of question: are extreme values out of place?)
    - [ ] Are structures from R/L more prone to be at the top or bottom? Similar kind of test 
    - [X] Account for stages in permutations in Sestan
----  
  
- [ ] Update **observations.md** file: Target date Friday 12


--- 
  
- [ ] Intregrate **code** into drake workflow
  + [ ] Check log2 consistency
  + [ ] Erase non-parametric tests if not used (superseeded)

--- 
--- 

- [X] **Stats** for SNP count and frequency (not modified)
- [X] **Stats** for cell type enrichment (not modified)
- [X] Mean expression DEG from Cell Atlas data (modified with log2)
  + [X] Stats (pairwise paired t tests) + heatmaps plot
  
- [X] **PCA** for the whole (414 structures) dataset.
  + Comparison to PCAs for 414 Akey and 414 AkeyPey  
  


