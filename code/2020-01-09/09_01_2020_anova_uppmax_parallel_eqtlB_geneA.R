library(data.table)
library(parallel)

effect_eqtl_gene <- function(res, pheno=phenotype, geno=genotype){
  # function that calculates the effect of a qtl on a gene - using anova
  # outputs a string with anova.pval__anova.r2
  # res - data.table with 2 columns; 
  # first column: gene 
  # second column: eqtl
  lmp <- function (modelobject) {
    # function to get the p-value out of a lm
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
}

start_time <- Sys.time()

cl = makeCluster(detectCores(), type="FORK")
path <- "/home/carolpb/DegreeProject/"
phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
effectsA_B.sepA_B <- fread(paste0(path,"results/2020-01-07/infoA_B.gz"))

eqtls.B <- na.omit(unique(effectsA_B.sepA_B$eqtl.B))
genes.A <- unique(effectsA_B.sepA_B$geneA)
res.tot.eqtlB_geneA <- data.table(expand.grid(gene=genes.A, eqtl=eqtls.B))#, anv.res=NA))
res.eqtlB_geneA <- res.tot.eqtlB_geneA[1:1000,]

clusterEvalQ(cl, c(library(data.table), library(parallel)))
clusterExport(cl, c("path", "phenotype", "genotype", "effectsA_B.sepA_B", "effect_eqtl_gene", "res.eqtlB_geneA"))

res.eqtlA_geneB$anv.res <- parApply(cl=cl,res.eqtlA_geneB,1,effect_eqtl_gene, phenotype, genotype)
stopCluster(cl)


fwrite(res.eqtlA_geneB, paste0(path, "results/2020-01-09/anova_eqtlB_geneA.gz"))
end_time <- Sys.time()
end_time - start_time


