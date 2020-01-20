# Functions I've created for my analysis + some useful ones found on the internet


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
  # merge the results of the anova with the main table
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
  
  geno_1 <- sum(ifelse(genotype[,..snp] == 1, 1,0))
  geno_neg1 <- sum(ifelse(genotype[,..snp] == -1, 1,0))
  
  return(paste(geno_1, geno_neg1, sep="__"))
}

flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations
  
  library("dplyr")
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
}

expand.grid.faster <- function(seq1,seq2) {
  # faster alternative to expand.grid
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}

create_ini_table <- function(eqtl.table, genesB, var.exp.lim){
  # requires data.table, tidyr
  # takes a talbe with 4 columns: gene, eqtl, cis and var.exp and creates a table with all the possible combinations
  # of geneA-eqtlA with geneB-eqtlB
  # cis is a True/false - if the eqtl is in cis with the gene
  # genesB is the genes we want to compare the gene-eqtl pairs with
  # var.exp.lim is chosen by the user
  # For the first set of pairs, only the gene-eqtl that are in cis and where the variance explained is 
  # above the set limit will be kept
  
  # Keep only the gene-eqtl pairs where the eqtl is in cis with the gene and where the var.exp >0.1
  eqtl.tableA <- eqtl.table[cis==T & var.exp > var.exp.lim]
  
  # unite columns so they act as one block of information
  eqtl.tableA.unite <- eqtl.tableA %>% unite(infoA, gene, pmarker, cis, var.exp, sep = "__")
  
  # get table with all genes that are going to be tested with eqtlA (geneB migth have an eqtl or not)
  eqtl.tableB<- merge(data.table(genesB),eqtl.table, by.x="genesB", by.y="gene", all.x=T)
  
  # unite columns so they act as one block of information
  eqtl.tableB.unite <- eqtl.tableB %>% unite(infoB, genesB, pmarker, cis, var.exp, sep = "__")
  
  # get all combinations of geneA and geneB with corresponding eqtls
  eqtl.table.combineAB <- expand.grid.faster(eqtl.tableA.unite$infoA,eqtl.tableB.unite$infoB)
  eqtl.table.combineAB.dt <- data.table(eqtl.table.combineAB)
  setnames(eqtl.table.combineAB.dt, old=c("Var1", "Var2"), new=c("infoA", "infoB"))
  
  # separate info blocks into normal columns again
  eqtl.table.sepA <- eqtl.table.combineAB.dt %>% separate(infoA, c("geneA", "eqtl.A", "cis.A", "var.exp.A"), sep="__")
  eqtl.table.sepAB <- eqtl.table.sepA %>% separate(infoB, c("geneB", "eqtl.B", "cis.B", "var.exp.B"), sep="__")
}

find.effects_fun <- function(effects.table, snp.pval, snp.pval.nsign){
  # the input table must have two columns called "eqtlA_geneB.pval" and "eqtlB_geneA.pval", which contain the p-val of the anova of eqtlX on geneY
  # This function does not take into account if the eqtl is in cis with the gene, so the input table must have been subseted to only contain eqtl-gene pairs that are in cis
  
  # adds two columns to the input table: A->B and B->A. 
  # These columns can have a value of 
  #   * NA- if with the chosen p-val cutoffs we can't say if X affects Y
  #   * TRUE - if X affects Y
  #   * FALSE - if X does not affect Y
  
  # if both genes have same eqtl, then it's the third case, C, where both are affected by the same. 
  # we can't say anything about the relationship between A and B. 
  # We just know that a third factor is affecting them
  
  effects.table[, ("A->B") := NA]; effects.table[as.numeric(eqtlA_geneB.pval) < snp.pval, ("A->B") := T];
  #effects.table[(as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) ), ("A->B") := F]
  effects.table[(as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B))), ("A->B") := F]
  
  effects.table[, ("B->A") := NA]; effects.table[as.numeric(eqtlB_geneA.pval) < snp.pval, ("B->A") := T];
  #effects.table[(as.numeric(eqtlB_geneA.pval) > as.numeric(snp.pval.nsign)), ("B->A") := F]
  effects.table[(as.numeric(eqtlB_geneA.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B))), ("B->A") := F]
  
  return(effects.table)
}

testparams <- function(params, tab){
  #function to get cases where T and F with different params
  # to run in parallel
  sign_p <- params[1]
  non_sign_p <- params[2]
  cor_p <- params[3]
  
  
  temp <- tab[cor.pval < cor_p & cis.A ==T & cis.B==T & geneA!=geneB]
  temp <- find.effects_fun(temp, sign_p, non_sign_p)
  
  temp_TF.1 <- temp[temp$`A->B`==T & temp$`B->A`==F, .(geneA, geneB, eqtl.A, eqtl.B, `A->B`, `B->A`)]
  temp_TF.2 <- rbind(temp_TF.1, temp[temp$`A->B`==F & temp$`B->A`==T, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
  temp_TF <- unique(temp_TF.2)
  
  temp_TF$sign.p <- sign_p
  temp_TF$non_sign_p <- non_sign_p
  temp_TF$cis_p <- cor_p
  
  return(temp_TF)
}