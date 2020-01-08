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


effectsA_B.sepA_B <- fread(paste0(path,"results/2020-01-07/infoA_B.gz"))



##### FUNCTIONS #####

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

merge_after_anova <- function(res.effect_eqtl_gene, gene.AB, eqtl.AB, effects.table){
  # requires library(tidyr) for the separate
  gene.AB <- toupper(gene.AB)
  eqtl.AB <- toupper(eqtl.AB)
  
  pval.name <- paste0("eqtl",eqtl.AB, "_", "gene", gene.AB, ".pval")
  r2.name <- paste0("eqtl",eqtl.AB, "_", "gene", gene.AB, ".r2")
  
  res.sep <- res.effect_eqtl_gene %>% separate(anv.res, c(pval.name, r2.name), sep="__")
  
  setnames(res.sep,"gene", paste0("gene", gene.AB))
  setnames(res.sep,"eqtl", paste0("eqtl.", eqtl.AB))
  
  effects.eqtlB <- merge(effects.table, res.sep, by=c(paste0("gene", gene.AB),paste0("eqtl.", eqtl.AB)), all.x=T)
  
  return(effects.eqtlB)
}

get_num_genos <- function(res, genotype){
  # get the number of each geno per snp
  # res is a data.table with one row - the snp ids
  snp <- res[1]
  
  geno_0 <- sum(ifelse(genotype[,..snp] == 0, 1,0)) # for one snp
  geno_1 <- sum(ifelse(genotype[,..snp] == 1, 1,0))
  geno_2 <- sum(ifelse(genotype[,..snp] == 2, 1,0))
  
  return(paste(geno_0,geno_1, geno_2, sep="__"))
}

#######
### Find effect of eqtls from geneA in expression of geneB
eqtls.A <- unique(effectsA_B.sepA_B$eqtl.A)
genes.B <- unique(effectsA_B.sepA_B$geneB)
res.tot.eqtlA_geneB <- data.table(expand.grid(gene=genes.B, eqtl=eqtls.A))#, anv.res=NA))
res.eqtlA_geneB <- res.tot.eqtlA_geneB[1:10,]

# run the anova function in parallel -- effect of eqtlA on geneB
system.time({
cl = makeCluster(detectCores() - 1, type="FORK")
res.eqtlA_geneB$anv.res <- parApply(cl=cl,res.eqtlA_geneB,1,effect_eqtl_gene, phenotype, genotype)
stopCluster(cl)
})

### Find effect of eqtls from geneB in expression of geneA
eqtls.B <- na.omit(unique(effectsA_B.sepA_B$eqtl.B))
genes.A <- unique(effectsA_B.sepA_B$geneA)
res.tot.eqtlB_geneA <- data.table(expand.grid(gene=genes.A, eqtl=eqtls.B))#, anv.res=NA))
res.eqtlB_geneA <- res.tot.eqtlB_geneA[1:10,]

# run the anova function in parallel -- effect of eqtlB on geneA
system.time({
  cl = makeCluster(detectCores() - 1, type="FORK")
  res.eqtlB_geneA$anv.res <- parApply(cl=cl,res.eqtlB_geneA,1,effect_eqtl_gene, phenotype, genotype)
  stopCluster(cl)
})


## merge anova results with the information table to create an effects table
# results from effect of eqtlA on geneB
effects_table.eqtlA_geneB <- merge_after_anova(res.eqtlA_geneB, "B", "A", effectsA_B.sepA_B)

# results from effect of eqtlB on geneA
effects_table.anova <- merge_after_anova(res.eqtlB_geneA, gene.AB="A", eqtl.AB="B", effects_table.eqtlA_geneB )




# Find number of each geno for each SNP
snp_ids <- colnames(genotype)[2:length(colnames(genotype))]
nSNPs <- length(snp_ids)

ngenos.dt <- data.table(snp_ids[1:1000])

system.time({
  cl = makeCluster(detectCores() - 1, type="FORK")
  ngenos.dt$nums <- parApply(cl=cl,ngenos.dt,1,get_num_genos, genotype)
  stopCluster(cl)
})

numgenos <- ngenos.dt %>% separate(nums, c("n0", "n1", "n2"), sep="__")
