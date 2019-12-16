### Anova of the effect of an eQTL on the gene it is affecting

library(data.table)
# path <- "/home/carolpb/DegreeProject/" # use with UPPMAX
path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" # use with own computer

phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
eqtl_results <- fread(paste0(path,"data/SI_Data_04_eQTL.csv"))


lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

eqtl_results$anova.pval <- NA
eqtl_results$anova.r2 <- NA

start_time <- Sys.time()
for (rownum in 1:nrow(eqtl_results)){
  gene <-  eqtl_results$gene[rownum]
  eqtl <- eqtl_results$pmarker[rownum]
  
  anv <- lm(unlist(phenotype[ , ..gene]) ~ factor(unlist(genotype[ ,..eqtl])))
  
  eqtl_results$anova.pval[rownum] <- lmp(anv)
  eqtl_results$anova.r2[rownum] <- summary(anv)$adj.r.squared
}
end_time <- Sys.time()
end_time - start_time

eqtl_results[cis=="VERDADEIRO"]$cis <- TRUE
eqtl_results[cis=="FALSO"]$cis <- FALSE
fwrite(eqtl_results, paste0(path,"results/2019-12-11/eqtl_results_effect_eqtl-gene.gz"), sep="\t")

test <- fread(paste0(path,"results/2019-12-11/eqtl_results_effect_eqtl-gene.gz"))
eqtl_results[order(-rank(anova.pval))]


