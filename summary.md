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


#### **[File with my functions](https://github.com/CarolinaPB/DegreeProject/blob/master/code/myfunctions.R)**

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
load(paste0(path, "results/2020-01-10/effectstable.Rdata"))
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
> [Script](https://github.com/CarolinaPB/DegreeProject/blob/master/code/2020-01-20/20_01_2020_testparams.R)

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

## 2020-01-24
Defined new parameters that will be used from now in the analysis. Did the causal analysis again.
> [Script](https://github.com/CarolinaPB/DegreeProject/blob/master/code/2020-01-24/24_01_2020_changingparams.R)


Changed parameters being used based on the previous plots

parameter | exp | value |
--|--|--|
var.exp.lim| | 0.1
nSNPs| | 42052
nGenes| | 5720
snp.pval| | 0.01 *
snp.pval.sign| | 1e-5 |
corr.pval | 0.05/choose(nGenes,2) | 3.056919e-09 *  

\* Values that are different than the original parameters



Number of times each gene is the causal one or is on the receiving end

Subset:
```
   gene     causal  receiving
YLR270W     897             5
YLR264W     834             6
YLR281C     752             0
YNL088W     742             3
YLR260W     739            10
```

| ![num gene pairs](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/numtimes_genepair_newparams.png) |
|:--:|
| *Number of unique gene pairs for the cases where A->B and not B->A when using the new parameters* |

With these new parameters, many more gene pairs are found. The number of pairs that appear once increases to more than double of what we had before, the number of pairs that appear twice or three times also increases and now we have gene pairs that appear four times

| ![new params blob](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/TFgenes_blob_newparams.png) |
|:--:|
| *Causality network for the cases where A->B and not B->A when using the new params* |

## 2020-01-24 -- 2020-01-28
Since all genes are connected in my network and there are no small causal clusters, I looked for ways to group the genes. I ended up using the link community method (Ahn et al., 2010) to find sets of genes that are more highly connected with each other than with the rest of the network. This method was applied throught the linkcomm R package (Kalinka & Tomancak 2011). Using this method, nodes (genes) may be present in more than one cluster.

> [Get link communities](https://github.com/CarolinaPB/DegreeProject/blob/master/code/2020-01-27/27_01_2020_netcomm.R) and [cluster analysis](https://github.com/CarolinaPB/DegreeProject/blob/master/code/2020-01-29/29_01_2020_clusteranalysis.Rmd)

| ![linkcomm summary](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/linkcommsummary.png) |
|:--:|
| *summary of the linkcomm analysis. 3040 communities were found* |

| ![community membership](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/communitymembership.png) |
|:--:|
| *community membership for nodes that belong to the most communities. Colours indicate community-specific membership. Numbers on top of the matrix are the number of the clusters; Numbers on the right are the number of clusters those genes belong to; numbers on the bottom are the number of genes in each cluster* |

The next figure shows two clusters that share genes. The fraction of the total number of edges that a node has in each community is depicted using a pie chart

| ![example cluster](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/examplecluster1and7.png) |
|:--:|
| *Cluster 1 and cluster 7 share 3 genes* |

There are several measures to describe the centrality of a cluster, such as modularity (high modular clusters have more links between the cluster's nodes than going out of the cluster) and connectivity (there are more links going out of the clusters than within)
Networks with high modularity have dense connections between the nodes within modules but sparse connections between nodes in different modules.
Modularity is often used in optimization methods for detecting community structure in networks. [wikipedia](https://en.wikipedia.org/wiki/Modularity_(networks))


| ![top10modularity](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/top10modularity.png) |
|:--:|
| *clusters with top 10 modularity - more connections within the cluster than outside the cluster* |


| ![top10modularityclusters](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/top10modularclusters.png) |
|:--:|
| *10 clusters that have the highest modularity* |

| ![modularity all clusters](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/clustermodularity.png) |
|:--:|
| *modularity for all clusters* |

Only a few of the clusters have very high modularity

| ![top10connectivity](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/top10connectivity.png) |
|:--:|
| *clusters with top 10 connectivity (the inverse of modularity)* |


| ![top10connectedclusters](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/top10connectedclusters.png) |
|:--:|
| *10 clusters that have the highest connectivity* |

There's still only one cluster with no subclusters

Comparing with the previous network plot, we can see that the clusters with high connectivity are all interconnected and have a lot of links going out of them. The clusters in the highest modularity network plot are mostly unconnected.


## 2020-02-05 -- 2020-02-12

In order to find the GO terms associated with my genes I used YeastMine (Balakrishnan et al., 2012) (https://yeastmine.yeastgenome.org/yeastmine/begin.do). Since I was not being able to do it in R, using the YeastMine API, I used python to run my queries. To be able to run the queries with python, first I needed to create and account and request an API key. Since you can generate python code from the website, I used it as a guide and added/ removed parameters to get what I needed.
From YeastMine I got the GO code and term for my genes, as well as if it belongs to a biological process or other, and the evidence code for the annotation.


> [GO analysis](https://github.com/CarolinaPB/DegreeProject/blob/master/code/2020-02-07/19_02_2020_GO_analysis.Rmd)

One gene might be associated with more than one GO term


| ![barplot top links GO](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/numlinksGOterm_genesA_top10.png) |
|:--:|
| *Causal genes GO terms that have the most links going out from them* |

| ![barplot top links GO](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/piechart_top10_GOterms_genesA.png) |
|:--:|
| *The top 10 most represented GO term categories in the causal genes* |

According to [EBI](https://www.ebi.ac.uk/QuickGO/term/GO:0050789), the GO term "biological process" is "any process that modulates the frequency, rate or extent of a biological process" it's usually used when the actual function of the gene is not known.  

High values might be because there are many genes associated with a certain term of because the genes that are associated with that term have many "arrows" going out (they affect many other genes)


Before looking at enrichment, we hypothesised that many causal genes would be associaded with transcription/transcription factors.



| ![boxplot transc factor](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/linkstranscriptionboxplot.png) | ![boxplot transc factor outliers](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/linkstranscriptionboxplot_outliers.png)
|:--:|:--:


How many links go from genes for with the GO term includes "transcription factor". Left - without outliers plotted

| ![boxplot transc factor/regulator](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/linkstranscriptionregulationboxplot.png) | ![boxplot transc factor/regulator outliers](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/linkstranscriptionregulationboxplot_outliers.png)
|:--:|:--:

How many links go from genes for with the GO term includes "transcription factor" or "transcription" and "regulation". Left - without outliers plotted

We were expecting that the boxplot including transcription factors/regulators would have a higher mean value than the other

# GO Enrichment
## 2020-02-13 -- 2020-02-14   

> [enrichment analysis](https://github.com/CarolinaPB/DegreeProject/blob/master/code/2020-02-13/13_02_2020_enrichmentanalysis_hypergeo.R)

I used GOstats, an R package (bioconductor), to test GO terms for over representation. I used both a classical hypergeometric test and a conditional hypergeometric test, which uses the relationships among GO terms to decorrelate the results

First I needed to define a few parameters:
* **universe** - all the genes in the dataset (can be involved in the causality or not) (num genes = 5720)
* **interesting genes** - causal genes (n=2658) or affected genes (n=2478)

Falcon & Gentleman (2007) the universe can be reduced by not using the genes that are not being expressed (in this case I would say not involved in the causality). Taking this into account, it would be interesting to perform the hypergeometric test using only the genes involved in the causality (genes that affect the expression of other genes and genes that are affected) as universe. Falcon & Gentleman (2007) also suggest removing genes that do not map to any GO term
I'm performing the hypergeometric test twice, once for the causal genes and once for the affected genes to see if there's a different enrichment in both groups. It would be expected that the causal group would be enriched for genes involved in regulation.

From Falcon & Gentleman (2007)  
"In the hypergeometric model, each term is treated as an independent classification. Each gene is classified according to whether or not it has been selected and whether or not it is annotated at a particular term. A hypergeometric probability is computed to assess whether the number of selected genes associated with the term is larger than expected."

Performed new hypergeometric test with
* **universe** - all the genes involved in the causality (num genes = 2861)
* **interesting genes** - causal genes or affected genes

I checked if the resulting enrichment table was the same for both universes tested and it was.
I will continue by using the results from the second test, where the universe was comprised of genes involved in the causality.

### Hypergeo results

The following overviews show the number of significant terms found
| ![hypergeo genesA summary](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/hypergeosummary_genesA.png) |
|:--:|
| *overview of hypergeometric test performed on the causal genes* |


From the 4425 "biological process" GO terms tested, 114 were found to be overenriched (at p-value 0.05)

| ![hypergeo genesB summary](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/hypergeosummary_genesB.png) |
|:--:|
| *overview of hypergeometric test performed on the affected genes* |

From the 4356 "biological process" GO terms tested, 114 were found to be overenriched (at p-value 0.05)

### Conditional Hypergeo results

| ![cond hypergeo genesA summary](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/condhypergeosummary_genesA.png) |
|:--:|
| *overview of conditional hypergeometric test performed on the causal genes* |


| ![cond hypergeo genesB summary](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/condhypergeosummary_genesB.png) |
|:--:|
| *overview of conditional hypergeometric test performed on the affected genes* |

set of genes | Hypergeo | Conditional hypergeo
-|-|-
causal | 114 | 51
affected | 114 | 67

To see the enrichment results for the causal genes (conditional hypergeo test) see this [table](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/hgCondA_htmlreport.md) (contains links to information on each GO term). For the affected genes' enrichment check this [table](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/hgCondB_htmlreport.md)

| ![graphA.1 overview](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/termgraphA_graph1_overview.png) | ![graphA.1 sub](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/termgraphA_graph1_subsection.png) | ![graphcondA.1](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/termgraphcondA_graph1.png)
|:--:|:--:|:--:|

Two first figures: one of the GO term graphs from the hypergeometric test for the causal genes.
Figure on the right: one of the graphs resulting from the conditional hypergeometric test.
The area delimited in blue is amplified on the center figure. These are terms that have as a parent term "RNA metabolic process". On the right, the "simplified" DAG where the parent term is "RNA metabolic process" after the conditional hypergeometric test.

(analysis inspired by (Hahne, Huber, Gentleman, & Falcon, 2008))  
Files where you can see a all the subgraphs of the gene ontology directed acyclic graph for the enrichment of the [causal genes](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/termgraph_A_bioproc_cond.pdf) and for the [affected genes](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/termgraph_B_bioproc_cond.pdf). The subgraphs consist of nodes that are connected in the DAG, and where all nodes are significant (according to the hypergeo test).   
The arrows in the plots point from a more specific term (child) to a more general term (parent).  
The conditional hypergeo test is trying to answer the question "Is there additional evidence to mark the a certain parent term significant beyond that provided by its significant child?"
As said before, the conditional hypergeo test uses the relationship between the GO terms to adjust the enrichment results. The method implemented by Hahne, Huber, Gentleman, & Falcon (2008) conditions on all child terms which are themselves significant at a specified p-value cutoff.


#### GO enrichment comparison for causal and affected genes
See the high quality image at [heatmap](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/heatmap_enrichmentpvals.pdf) (download to be able to zoom in and read the name of the GO terms better)

| ![heatmap overview](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/pheatmap_overview_enrichmentpval_negativelog10.png) |
|:--:|
| *Heatmap of GO enrichment p-values (-log10)* |

The left column represents the terms associated with the causal genes and the right column with the affected genes. There's a clear separation of terms between causal and affected genes. We can see that there's only one shared GO term between causal and affected genes.  
In the next two figures we can see that 25 out of 51 terms associated with the causal genes are involved in regulation of some process. The fact that the causal genes are enriched for terms involved in regulation and that the affected genes are not makes it seem like our way to infer causality is working or is going in the right direction.  

| ![heatmap sub1](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/pheatmap_enrichmentnegativelog10pval_sub1.png) |
|:--:|
| *First section of the heatmap of GO enrichment p-values (-log10)* |

| ![heatmap sub2](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/pheatmap_enrichmentnegativelog10pval_sub2.png) |
|:--:|
| *Second section of the heatmap of GO enrichment p-values (-log10)* |



### Enrichment with different genes of interest
Counted the genes going in or out of a gene and plotted the distribution:

| ![scatter linknum](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/scatter_linknum.png) |
|:--:|

| ![scatter linknum zoom](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/scatter_linknum_zoom.png) |
|:--:|

We can see that if the genes have a high number of "links" going out, they don't have many going "in" and that if they have less than around 40 links out, they have a high number of links going in. This shows that a clear separation between genes that affect other genes and genes that are being affected.


| ![heatmap >40 links causal](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/heatmap_enrichmentpvals_geneslinks_c40.png) |
|:--:|
| *Heatmap of GO enrichment p-values (-log10) for the causal genes that affect at least 40 other genes and for the affected gnes that don't affect any other gene* |

Click [here](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/heatmap_enrichmentpvals_geneslinks_c40.pdf) for the pdf - download to be able to zoom.

| ![heatmap >40 links causal](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/heatmap_enrichment_geneslinks_c40_subcausal.png) |
|:--:|
| *Section corresponding to the causal genes* |

| ![heatmap >40 links causal](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/heatmap_enrichment_geneslinks_c40_subaffected.png) |
|:--:|
| *Section corresponding to the affected genes* |


### Plot with positions of affected gene vs positions of causal genes
| ![genepos vs genepos](https://github.com/CarolinaPB/DegreeProject/blob/master/results/results_figures/images/affectedgenes_vs_causalgenes_position.png) |
|:--:|
| *Affected gene position vs causal gene position* |

Several of the vertical "bands" match with the ones reported in Albert et al (2018) - Figure 3.  
In my figure I'm calling "bands" to cases where the same gene (x axis) is affecting several genes (y axis).  
I seem to have bands in the same chromosomes where Albert et al (2018) report eQTL hotspots - for example, the "bands" on chr VII, XII, XIV
