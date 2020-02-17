library(data.table)
library(tidyr)

path <- "/home/carolpb/DegreeProject/" # use with UPPMAX

phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
infoA_B <- fread(paste0(path, "results/2020-01-07/infoA_B.gz"))

effect_eqtl_gene <- function(res, pheno=phenotype, geno=genotype){
  
  lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }
  
  gene <- res[1]
  eqtl <- res[2]
  
  anv <- lm(unlist(pheno[ , ..gene]) ~ unlist(geno[ ,..eqtl]))
  anv.res <- paste(lmp(anv), summary(anv)$adj.r.squared, sep="__")
  
  return(anv.res)
  # return(gene)
}

res$anv.res <- apply(res,1,effect_eqtl_gene, phenotype, genotype)

effect_eqtlA_geneB <- res %>% separate(anv.res, c("eqtlA_geneB.pval", "eqtlA_geneB.r2"), sep="__")
