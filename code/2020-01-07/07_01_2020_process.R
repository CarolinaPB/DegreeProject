##### LOAD FILES #####
library(data.table)
library(tidyr)
library(parallel)

# path <- "/home/carolpb/DegreeProject/" # use with UPPMAX
path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" # use with own computer

phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
#eqtl_results <- fread(paste0(path,"results/2019-12-11/eqtl_results_effect_eqtl-gene.gz"))
eqtl_results <- fread(paste0(path,"data/SI_Data_04_eQTL.csv"))
eqtl_results[cis=="VERDADEIRO"]$cis <- T
eqtl_results[cis=="FALSO"]$cis <- F

##### PARAMETERS #####
var.exp.lim <- 0.1

#####
#subset eqtl_results table
eqtl_results.sub <- eqtl_results[,.(gene, pmarker, cis, var.exp)]

# the first group of information in the table we are creating will correspond to the genes that have the eqtl in cis and that have var.exp > lim
# Keep only the gene-eqtl pairs where the eqtl is in cis with the gene and where the var.exp >0.1
eqtl_results.sub.cis <- eqtl_results.sub[cis==T & var.exp > var.exp.lim]

# unite columns so they act as one block of information
eqtl_results.sub.uniteA <- eqtl_results.sub.cis %>% unite(infoA, gene, pmarker, cis, var.exp, sep = "__")

# the second group of information in the table corresponds to all the genes (with or without eqtl)
genesB <- colnames(phenotype[,2:ncol(phenotype)])

# get table with all genes that are going to be tested with eqtlA (geneB migth have an eqtl or not)
eqtl_results.sub.B<- merge(data.table(genesB),eqtl_results.sub, by.x="genesB", by.y="gene", all.x=T)

# unite columns so they act as one block of information
eqtl_results.sub.uniteB <- eqtl_results.sub.B %>% unite(infoB, genesB, pmarker, cis, var.exp, sep = "__")

# get all combinations of geneA and geneB with corresponding eqtls
effectsA_B.temp <- expand.grid(infoA =eqtl_results.sub.uniteA$infoA,infoB=eqtl_results.sub.uniteB$infoB)
effectsA_B.temp.dt <- data.table(effectsA_B.temp)

# separate info blocks into normal columns again
effectsA_B.sepA <- effectsA_B.temp.dt %>% separate(infoA, c("geneA", "eqtl.A", "cis.A", "var.exp.A"), sep="__")
effectsA_B.sepA_B <- effectsA_B.sepA %>% separate(infoB, c("geneB", "eqtl.B", "cis.B", "var.exp.B"), sep="__")

# Final table with geneA and geneB and their corresponding eqtls (geneB might not have an eqtl)
effectsA_B.sepA_B

fwrite(effectsA_B.sepA_B,"results/2020-01-07/infoA_B.gz")



### Find effect of eqtls from geneA in expression of geneB
eqtls <- unique(effectsA_B.sepA_B$eqtl.A)
genes <- unique(effectsA_B.sepA_B$geneB)
res.tot <- data.table(expand.grid(gene=genes, eqtl=eqtls, anv.res=NA))
res <- res.tot[1:2000,]

lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


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

# res$anv.res <- apply(res,1,effect_eqtl_gene)

system.time({
cl = makeCluster(detectCores() - 1, type="FORK")
res$anv.res <- parApply(cl=cl,res,1,effect_eqtl_gene, phenotype, genotype)
stopCluster(cl)
})



effect_eqtlA_geneB <- res %>% separate(anv.res, c("eqtlA_geneB.pval", "eqtlA_geneB.r2"), sep="__")







# for (rownum in 1:nrow(res)){
#   gene <- res$gene[rownum]
#   eqtl <- res$eqtl[rownum]
#   
#   anv <- lm(unlist(phenotype[ , ..gene]) ~ unlist(genotype[ ,..eqtl]))
#   res$anova.pval[rownum] <- lmp(anv)
#   res$anova.r2[rownum] <- summary(anv)$adj.r.squared
# }



# TODO
# get effect of eqtlA in all unique geneB --> function









gene <-  eqtl_results$gene[rownum]
eqtl <- eqtl_results$pmarker[rownum]

anv <- lm(unlist(phenotype[ , ..gene]) ~ unlist(genotype[ ,..eqtl]))

eqtl_results$anova.pval[rownum] <- lmp(anv)
eqtl_results$anova.r2[rownum] <- summary(anv)$adj.r.squared

