#   ____________________________________________________________________________
#   Load stuff                                                              ####
library(data.table)
path <- "/home/carolpb/DegreeProject/" # use with UPPMAX
#path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" # use with own computer

phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
eqtl_results <- fread(paste0(path,"data/SI_Data_04_eQTL.csv"))


#   ____________________________________________________________________________
#   Functions                                                               ####
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


#   ____________________________________________________________________________
#   Commands                                                                ####
block_num = commandArgs(trailingOnly=TRUE)
message(block_num)

eqtls.list <- unique(eqtl_results$pmarker)
allgenes <- colnames(phenotype[,2:ncol(phenotype)])

chuncksize <- 1000
block <- as.numeric(block_num)
start <- (block*chuncksize) +1
end <- (start + chuncksize) -1
end <- min(end, length(allgenes))

allgenes <- allgenes[start:end]

gene_eqtl.anova <- data.table(matrix(ncol = 4, nrow = length(eqtls.list)*length(allgenes)))
colnames(gene_eqtl.anova) <- c("gene", "eqtl", "pval", "r2")

start_time <- Sys.time()
rown <- 1
for (eqtl in eqtls.list){
  for (gene in allgenes){
    gene_eqtl.anova$gene[rown] <- gene
    gene_eqtl.anova$eqtl[rown] <- eqtl
    
    anv <- lm(unlist(phenotype[ , ..gene]) ~ factor(unlist(genotype[ ,..eqtl])))
    gene_eqtl.anova$pval[rown] <- lmp(anv)
    gene_eqtl.anova$r2[rown] <- summary(anv)$adj.r.squared
    
    rown <- rown+1
    print(rown)
  }
}
end_time <- Sys.time()
end_time - start_time

fwrite(gene_eqtl.anova, paste0("/home/carolpb/DegreeProject/results/2019_12_09/09_12_2019_anova_", start, "_", end), sep="\t")
