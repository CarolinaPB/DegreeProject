# Summary of what I've been doing

Started by playing around with the data.  

My starting dataset:
```r
phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
eqtl_results <- fread(paste0(path,"data/SI_Data_04_eQTL.csv"))
```

Data from **Albert FW, Bloom JS, Siegel J, Day L, Kruglyak L. 2018. Genetics of trans-regulatory variation in gene expression. eLife 7: 1–39.**  [Source data](https://elifesciences.org/articles/35471/figures#supp1)  

I started with :
* **phenotype matrix** - contains the gene expression data (expression levels in units of log2(TPM) for all genes and segregants)
* **genotype matrix** - contains the genotype information (genotypes at 42,052 markers for all segregants. BY (i.e. reference) alleles are denoted by ‘−1’. RM alleles are denoted by ‘1’.)
* **eqtl_results** - Genes with a local eQTL and significant Allele-specific expression (ASE), and discordant direction of effect. (1) Positive values indicate higher expression in RM compared to BY. (2) Shown is the less sig- nificant p-value from the two ASE datasets. (3) The table shows only genes where both ASE datasets agreed in the direction of effect. Shown is the average effect.


##### **[File with my functions](https://github.com/CarolinaPB/DegreeProject/blob/master/code/myfunctions.R)**

## 2020-01-07 -- 2020-01-09
> [Summary script](https://github.com/CarolinaPB/DegreeProject/blob/master/code/2020-01-07/07_01_2020_process.R)  

Defined first set of parameters to try:

parameter | exp | value |
--|--|--|
var.exp.lim| | 0.1
nSNPs| | 42052
nGenes| | 5720
snp.pval| 0.05/(nGenes * nSNPs) | 2.078678e-10
snp.pval.sign| | 1e-5 |
corr.pval | 0.05/(nGenes*nGenes) | 1.528192e-09


From **eqtl_results** we have pairs of gene-eqtl and the information if the eqtl is in cis with the gene.
Defined **genesA** as a set of genes that have an eqtl and where the eqtl is in cis with the gene. The variance explained should also be over the var.exp.lim. **genesB** are all the genes that have expression data in **phenotype** (they can have an associated eqtl or not).
**eqtlA** is an eqtl in cis with geneA and **eqtlB** is an eqtl paired with geneB (not necessarily in cis)


Created table with all the possible combinations of pairings of geneA-eqtlA and geneB-eqtlB

```r
effectsA_B.sepA_B <- fread(paste0(path,"results/2020-01-07/infoA_B.gz"))
```

```
     geneA          eqtl.A cis.A   var.exp.A   geneB
1: YAL062W  chrI:33293_A/T  TRUE   0.5747081   Q0140
2: YAL060W  chrI:35818_A/G  TRUE   0.3220796   Q0140
3: YAL056W  chrI:39162_C/G  TRUE   0.6566113   Q0140
4: YAL054C  chrI:44384_A/G  TRUE   0.1493313   Q0140
5: YAL049C  chrI:52951_G/T  TRUE   0.7291269   Q0140

           eqtl.B   cis.B   var.exp.B    
chrIII:204897_G/C   FALSE   0.03919222
chrIII:204897_G/C   FALSE   0.03919222
chrIII:204897_G/C   FALSE   0.03919222
chrIII:204897_G/C   FALSE   0.03919222
chrIII:204897_G/C   FALSE   0.03919222
```

Used ANOVA to find the effect of a gene's eqtl on the other gene - effect of eqtlA on geneB and eqtlB on geneA

> [Summary script](https://github.com/CarolinaPB/DegreeProject/blob/master/code/2020-01-07/07_01_2020_process.R) and [do anova directory](https://github.com/CarolinaPB/DegreeProject/tree/master/code/2020-01-09).

```r
effects_table.anova <- fread(paste0(path, "results/2020-01-09/09_01_2020_anovatable.gz"))
```

Got correlation between all the genes and added it to the table with the anova results
```r
effects_table.cor <- fread(paste0(path, "results/2020-01-10/effectstable.gz"))
```

## Inferring causality

Assumptions:
* geneA is in cis with eqtlA
* geneB is in cis with eqtlB
* var.explained for geneA must be > var.exp.lim
* correlation pval is < corr.pval
Inferred if gene A is affecting geneB or if geneB is affecting geneA.
* geneA != geneB

There are two categories and several end results:  
Categories:
* A affects B: A->B
* B affects A: B->A

End results:
* **A->B = T and B->A = F** or **A->B = F and B->A = T** --> this is the case we are mostly interested in. It means we can say that a gene affects the other, but it's not affected by it.
* **A->B = T and B->A = NA** or **A->B = NA and B->A = T** --> we can say that a gene affects the other, but we can't say if the second gene affects the first
* **A->B = NA and B->A = NA** --> we can't say anything about causality
* **A->B = F and B->A = T** or **A->B = F and B->A = T** --> neither gene affects the other
* **A->B = T and B->A = T** or **A->B = T and B->A = T**

How it works:  
* **A->B = T** if anova p-value for the effect of eqtlA on geneB is < snp.pval
* **A->B = F** if the anova p-value of the effect of eqtlA on geneB is > snp.pval.nsign and geneA and geneB have different eqtls
* **B->A = T** if anova p-value for the effect of eqtlB on geneA is < snp.pval
* **B->A = F** if the anova p-value of the effect of eqtlB on geneA is > snp.pval.nsign and geneA and geneB have different eqtls

Table with results for A->B and B->A for all cases
```r
find.effects <- fread(paste0(respath, "2020-01-27/findeffects_all.gz"))
```

```
  geneA     geneB               eqtl.A            eqtl.B cis.A cis.B var.exp.A  var.exp.B
YAL003W   YAL008W      chrI:133174_G/A   chrI:136961_T/C  TRUE  TRUE 0.1190554 0.04176901
YAL003W   YAL009W      chrI:133174_G/A   chrI:132723_G/A  TRUE  TRUE 0.1190554 0.06344857
YAL003W   YAL010C      chrI:133174_G/A   chrI:134219_C/T  TRUE  TRUE 0.1190554 0.03777914
YAL003W   YAL013W      chrI:133174_G/A   chrI:131539_C/T  TRUE  TRUE 0.1190554 0.01516904
YAL003W   YAL022C      chrI:133174_G/A   chrI:114628_G/T  TRUE  TRUE 0.1190554 0.04732674

eqtlA_geneB.pval    eqtlA_geneB.r2  eqtlB_geneA.pval    eqtlB_geneA.r2         cor     cor.pval  A->B  B->A
2.857498e-08        0.029090909         6.810904e-15    0.057379498     -0.5419907 0.000000e+00    NA  TRUE
6.393194e-11        0.040463429         8.727858e-15    0.056923573     -0.4344359 0.000000e+00  TRUE  TRUE
4.723889e-06        0.019565417         5.127615e-14    0.053663535     -0.4371749 0.000000e+00    NA  TRUE
6.275525e-03        0.006389033         5.829682e-15    0.057665406     -0.4007383 0.000000e+00 FALSE  TRUE
4.631097e-07        0.023892823         9.575469e-08    0.026833917     -0.2238410 5.864198e-13    NA    NA

```


Table with with the gene-eqtl pairs where A->B =T and B->A=F

```r
find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))
```
```
geneA   geneB              eqtl.A             eqtl.B    A->B  B->A
YAL003W YAL033W   chrI:133174_G/A   chrI:84112_T/C      TRUE FALSE
YAL003W YAR050W   chrI:133174_G/A   chrI:202732_G/A     TRUE FALSE
YAL003W YDR059C   chrI:133174_G/A   chrIV:567417_C/T    TRUE FALSE
YAL003W YLL024C   chrI:133174_G/A   chrXII:86595_T/C    TRUE FALSE
YAL032C YAL049C   chrI:84112_T/C    chrI:52951_G/T      TRUE FALSE
```
#### Results using the parameters above
| ![initial blob](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/TFgenes_blob_firstparams.png) |
|:--:|
| *Causality network for the cases where A->B and not B->A* |

There are no sub clusters, only one unique network where all genes are connected

| ![num gene pairs](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/numtimes_genepair_firstparams.png) |
|:--:|
| *Number of unique gene pairs for the cases where A->B and not B->A* |

Three plots basically showing the same

| ![num gene-eqtl pairs](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/numtimes_geneeqtl_pairs_firstparams.png) |
|:--:|
| *Number of unique gene-eqtl pairs for the cases where A->B and not B->A* |

| ![num gene-eqtl pairs histogram](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/numtimes_geneeqtl_pairs_firstparams_hist.png) |
|:--:|
| *Number of unique gene-eqtl pairs for the cases where A->B and not B->A* |

| ![num gene-eqtl pairs histogram](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/TNAgenes_blob_firstparams.png) |
|:--:|
| *Causality network for the cases where A->B and we can't say if B->A or not* |

It seems like there's a large cluster of genes that are affecting each other, as well as smaller isolated clusters


## 2020-01-20
Test different parameter values combinations to see if I there is a certain combination that gives optimal results in the causality inference

Groups of parameters to test:  
```r
sign_p <- c(1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
non_sign_p <- c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
cor_p <- c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
```


Subset:
```
sign_p non_sign_p   cor_p
1e-17       1e-07   1e-09
1e-16       1e-07   1e-09
1e-15       1e-07   1e-09
1e-14       1e-07   1e-09
1e-13       1e-07   1e-09
```


| ![numgenepairs diff params](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/numgenepairs_bypval_testparams.png) |
|:--:|
| *Number of unique gene pairs found when using different cutoffs. Each pannel corresponds to a different correlation p-value cutff. X-axis is the -log(pvalue) (for the effect to be significant) and differenc colors represent different values for the non-significant p-value (for the effect to be considered non-significant)* |


| ![numgenepairs diff params](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/numgenepairs_bypval_testparams_nocorrpval.png) |
|:--:|
| *Number of unique gene pairs found when using different cutoffs. X-axis is the -log(pvalue) (for the effect to be significant) and differenc colors represent different values for the non-significant p-value (for the effect to be considered non-significant)* |
