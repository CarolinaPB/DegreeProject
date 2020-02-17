# Summary of what I've been doing

Started by playing around with the data.  

My starting dataset:
```{r}
phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
eqtl_results <- fread(paste0(path,"data/SI_Data_04_eQTL.csv"))
```

Data from **Albert FW, Bloom JS, Siegel J, Day L, Kruglyak L. 2018. Genetics of trans-regulatory variation in gene expression. eLife 7: 1–39.**  
I started with :
* **phenotype matrix** - contains the gene expression data (expression levels in units of log2(TPM) for all genes and segregants)
* **genotype matrix** - contains the genotype information (genotypes at 42,052 markers for all segregants. BY (i.e. reference) alleles are denoted by ‘−1’. RM alleles are denoted by ‘1’.)
* **eqtl_results** - Genes with a local eQTL and significant Allele-specific expression (ASE), and discordant direction of effect. (1) Positive values indicate higher expression in RM compared to BY. (2) Shown is the less sig- nificant p-value from the two ASE datasets. (3) The table shows only genes where both ASE datasets agreed in the direction of effect. Shown is the average effect.

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
