##### LOAD FILES #####
library(data.table)
library(tidyr)
library(parallel)
library(igraph)

# path <- "/home/carolpb/DegreeProject/" # use with UPPMAX
path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" # use with own computer
# path <- "/home/carolina/DegreeProject/" # use with diprotodon

phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
#eqtl_results <- fread(paste0(path,"results/2019-12-11/eqtl_results_effect_eqtl-gene.gz"))
eqtl_results <- fread(paste0(path,"data/SI_Data_04_eQTL.csv"))
eqtl_results[cis=="VERDADEIRO"]$cis <- T
eqtl_results[cis=="FALSO"]$cis <- F

##### PARAMETERS #####
var.exp.lim <- 0.1

nSNPs <- length(colnames(genotype))-1
nGenes <- length(colnames(phenotype))-1

snp.pval <- 0.05/(as.numeric(nGenes) * as.numeric(nSNPs))
snp.pval.nsign <- as.numeric(1e-5)

corr.pval <-  0.05/(nGenes*nGenes)

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

#######
### Find effect of eqtls from geneA in expression of geneB
eqtls.A <- unique(effectsA_B.sepA_B$eqtl.A)
genes.B <- unique(effectsA_B.sepA_B$geneB)
res.tot.eqtlA_geneB <- data.table(expand.grid(gene=genes.B, eqtl=eqtls.A))#, anv.res=NA))
res.eqtlA_geneB <- res.tot.eqtlA_geneB[1:1000,]

# run the anova function in parallel -- effect of eqtlA on geneB

start_time <- Sys.time()
  cl = makeCluster(detectCores() - 1, type="FORK")
  res.eqtlA_geneB$anv.res <- parApply(cl=cl,res.eqtlA_geneB,1,effect_eqtl_gene, phenotype, genotype)
  stopCluster(cl)
end_time <- Sys.time()
end_time - start_time


### Find effect of eqtls from geneB in expression of geneA
eqtls.B <- na.omit(unique(effectsA_B.sepA_B$eqtl.B))
genes.A <- unique(effectsA_B.sepA_B$geneA)
res.tot.eqtlB_geneA <- data.table(expand.grid(gene=genes.A, eqtl=eqtls.B))#, anv.res=NA))
res.eqtlB_geneA <- res.tot.eqtlB_geneA[1:1000,]

# run the anova function in parallel -- effect of eqtlB on geneA
start_time <- Sys.time()
  cl = makeCluster(detectCores() - 1, type="FORK")
  res.eqtlB_geneA$anv.res <- parApply(cl=cl,res.eqtlB_geneA,1,effect_eqtl_gene, phenotype, genotype)
  stopCluster(cl)
end_time <- Sys.time()
end_time - start_time


## merge anova results with the information table to create an effects table
# results from effect of eqtlA on geneB
effects_table.eqtlA_geneB <- merge_after_anova(res.eqtlA_geneB, "B", "A", effectsA_B.sepA_B)

# results from effect of eqtlB on geneA
effects_table.anova <- merge_after_anova(res.eqtlB_geneA, gene.AB="A", eqtl.AB="B", effects_table.eqtlA_geneB )
setcolorder(effects_table.anova, c("geneA", "geneB", "eqtl.A", "eqtl.B", "cis.A", "cis.B"))

fwrite(effects_table.anova, paste0(path, "results/2020-01-08/08_01_2020_anovatable_test.gz"), na = NA)
effects_table.anova <- fread(paste0(path, "results/2020-01-08/08_01_2020_anovatable_test.gz"))

# Find number of each geno for each SNP ######
snp_ids <- colnames(genotype)[2:length(colnames(genotype))]
nSNPs <- length(snp_ids)

ngenos.dt <- data.table(snp_ids)

system.time({
  cl = makeCluster(detectCores() - 1, type="FORK")
  ngenos.dt$nums <- parApply(cl=cl,ngenos.dt,1,get_num_genos, genotype)
  stopCluster(cl)
})

numgenos <- ngenos.dt %>% separate(nums, c("n1", "n-1"), sep="__")

fwrite(numgenos, paste0(path, "results/2020-01-08/08_01_2020_numgenos.gz"))

# all genos are represented by a large amount of samples ####



### correlation ####
cor_traits <- rcorr(as.matrix(phenotype[,2:ncol(phenotype)])) # to remove the sample names
save(cor_traits,file= "/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-01-08/full_cor.Rdata")

cor_traits.cor <- cor_traits$r
cor_traits.p <- cor_traits$P

cor_traits.cor[upper.tri(cor_traits.cor)] <- NA
cor_traits.p[upper.tri(cor_traits.p)] <- NA

my_cor_matr_flat <- flat_cor_mat(cor_traits.cor,cor_traits.p)
my_cor_matr_flat.dt <- data.table(my_cor_matr_flat)
my_cor_matr_flat.noNA <- na.omit(my_cor_matr_flat.dt)
colnames(my_cor_matr_flat.noNA) <- c("Trait1", "Trait2", "cor", "pval")

fwrite(my_cor_matr_flat.noNA, "/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-01-08/flat_cor_matr.gz")
#####
my_cor_matr_flat.noNA <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-01-08/flat_cor_matr.gz")


setnames(my_cor_matr_flat.noNA, "Trait1", "geneA")
setnames(my_cor_matr_flat.noNA, "Trait2", "geneB")
effects_table_cor.1 <- merge(effects_table.anova, my_cor_matr_flat.noNA, by=c("geneA", "geneB"), all.x = T)







######
# find.effects$`A->B` <- NA
# find.effects$`B->A` <- NA
# 
# find.effects[as.numeric(eqtlA_geneB.pval) < snp.pval]$`A->B` <-  TRUE
# find.effects[as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B))]$`A->B` <-  FALSE
# 
# find.effects[as.numeric(eqtlB_geneA.pval) < snp.pval]$`B->A` <-  TRUE
# find.effects[as.numeric(eqtlB_geneA.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.B) != as.character(eqtl.A))]$`B->A` <-  FALSE

find.effects_fun <- function(effects.table, dir){
  if (dir == "A_B"){
    dir.val="A->B"
  } else if (dir == "B_A"){
    dir.val="A->B"
  }
  effects.table[, (dir.val):=NA][as.numeric(eqtlA_geneB.pval) < snp.pval, (dir.val):=TRUE]
  effects.table[, (dir.val):=NA][as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B)), (dir.val):=FALSE]
  # effects.table[,(dir.val):=as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B))]
  # 
  # find.effects[as.numeric(eqtlA_geneB.pval) < snp.pval]$`A->B` <-  TRUE
  # find.effects[as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B))]$`A->B` <-  FALSE
}


# only keep the rows where there i
effects_table_nopvalNA <- effects_table.anova[!is.na(eqtlA_geneB.pval) | !is.na(eqtlB_geneA.pval)]

find.effects <- data.table(effects_table.anova) # use data.table to change the memory addresses (can be checked by tracemem(find.effects)==tracemem(effects_table.anova))


find.effects2[as.numeric(eqtlA_geneB.pval) < snp.pval]$`A->B` <-  TRUE

find.effects[, (dir.val):=NA][as.numeric(eqtlA_geneB.pval) < snp.pval, (dir.val):=TRUE]
find.effects[, (dir.val):=NA][as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B)), (dir.val):=FALSE]

find.effects[,(dir.val):=NULL]








effects <- find.effects

# one gene affects the other but we don't know the action of the other
case1 <- effects[(is.na(effects$`A->B`) & effects$`B->A` == TRUE) | (is.na(effects$`B->A`) & effects$`A->B` == TRUE)]
# one gene affects the other but is not affected by the other gene
case2 <- effects[(effects$`A->B` == FALSE & effects$`B->A` == TRUE) | (effects$`B->A` == FALSE & effects$`A->B` == TRUE)]
# we don't know if the genes affect each other
case3 <- effects[is.na(effects$`A->B`) & is.na(effects$`B->A`)]
# a gene affects the other but it's also affected by that other gene
case4 <- effects[(effects$`A->B` == TRUE & effects$`B->A` == TRUE) | (effects$`B->A` == TRUE & effects$`A->B` == TRUE)]
# no gene affects the other
case5 <- effects[(effects$`A->B` == FALSE & effects$`B->A` == FALSE) | (effects$`B->A` == FALSE & effects$`A->B` == FALSE)]

case1.A_to_B <- case1[case1$`A->B`==TRUE & is.na(case1$`B->A`)]
case1.A_to_B <- case1.A_to_B[,c("geneA","geneB", "A->B")]

case1.B_to_A <- case1[case1$`B->A`==TRUE & is.na(case1$`A->B`)]
case1.B_to_A <- case1.B_to_A[,c("geneB", "geneA", "B->A")]

if (nrow(case1.A_to_B) > 0 & nrow(case1.B_to_A) > 0){
  case1.links <- rbind(case1.A_to_B, case1.B_to_A, use.names=FALSE)
  
} else if (nrow(case1.A_to_B) == 0 & nrow(case1.B_to_A) > 0){
  case1.links <- case1.B_to_A
  
} else if (nrow(case1.A_to_B) > 0 & nrow(case1.B_to_A) == 0){
  case1.links <- case1.A_to_B
  
} else if (nrow(case1.A_to_B) == 0 & nrow(case1.B_to_A) == 0){
  case1.links <- case1.A_to_B
}

colnames(case1.links) <- c("from", "to", "type")
case1.links$case <- 1
case1.links <- unique(case1.links)

case2.A_to_B <- case2[case2$`A->B`==TRUE & case2$`B->A`==FALSE]
case2.A_to_B <- case2.A_to_B[,c("geneA","geneB", "A->B")]

case2.B_to_A <- case2[case2$`B->A`==TRUE & case2$`A->B`==FALSE]
case2.B_to_A <- case2.B_to_A[,c("geneB", "geneA", "B->A")]

if (nrow(case2.A_to_B) > 0 & nrow(case2.B_to_A) > 0){
  case2.links <- rbind(case2.A_to_B, case2.B_to_A, use.names=FALSE)
  
} else if (nrow(case2.A_to_B) == 0 & nrow(case2.B_to_A) > 0){
  case2.links <- case2.B_to_A
  
} else if (nrow(case2.A_to_B) > 0 & nrow(case2.B_to_A) == 0){
  case2.links <- case2.A_to_B
  
} else if (nrow(case2.A_to_B) == 0 & nrow(case2.B_to_A) == 0){
  case2.links <- case2.A_to_B
}

colnames(case2.links) <- c("from", "to", "type")
case2.links$case <- 2
case2.links <- unique(case2.links)

case3.links <- case3[,c("geneA", "geneB")]
case3.links$type <- NA
case3.links <- unique(case3.links)
colnames(case3.links) <- c("from", "to", "type")
case3.links$case <- 3

case4.A_to_B <- case4[case4$`A->B`==TRUE & case4$`A->B`==TRUE]
case4.A_to_B <- case4.A_to_B[,c("geneA","geneB", "A->B")]

case4.B_to_A <- case4[case4$`B->A`==TRUE & case4$`B->A`==TRUE]
case4.B_to_A <- case4.B_to_A[,c("geneB", "geneA", "B->A")]

if (nrow(case4.A_to_B) > 0 & nrow(case4.B_to_A) > 0){
  case4.links <- rbind(case4.A_to_B, case4.B_to_A, use.names=FALSE)
  
} else if (nrow(case4.A_to_B) == 0 & nrow(case4.B_to_A) > 0){
  case4.links <- case4.B_to_A
  
} else if (nrow(case4.A_to_B) > 0 & nrow(case4.B_to_A) == 0){
  case4.links <- case4.A_to_B
  
} else if (nrow(case4.A_to_B) == 0 & nrow(case4.B_to_A) == 0){
  case4.links <- case4.A_to_B
}

colnames(case4.links) <- c("from", "to", "type")
case4.links$case <- 4
case4.links <- unique(case4.links)


case5.A_to_B <- case5[case5$`A->B`==FALSE & case5$`A->B`==FALSE]
case5.A_to_B <- case5.A_to_B[,c("geneA","geneB", "A->B")]

case5.B_to_A <- case5[case5$`B->A`==FALSE & case5$`B->A`==FALSE]
case5.B_to_A <- case5.B_to_A[,c("geneB", "geneA", "B->A")]

if (nrow(case5.A_to_B) > 0 & nrow(case5.B_to_A) > 0){
  case5.links <- rbind(case5.A_to_B, case5.B_to_A, use.names=FALSE)
  
} else if (nrow(case5.A_to_B) == 0 & nrow(case5.B_to_A) > 0){
  case5.links <- case5.B_to_A
  
} else if (nrow(case5.A_to_B) > 0 & nrow(case5.B_to_A) == 0){
  case5.links <- case5.A_to_B
  
} else if (nrow(case5.A_to_B) == 0 & nrow(case5.B_to_A) == 0){
  case5.links <- case5.A_to_B
}

colnames(case5.links) <- c("from", "to", "type")
case5.links$case <- 5
case5.links <- unique(case5.links)



allcases.links <- do.call("rbind", list(case1.links, case2.links, case3.links, case4.links, case5.links))




#### Case1 plot
allcases.graph <- graph_from_data_frame(allcases.links[allcases.links$case == 1],directed=TRUE)

E(allcases.graph)$color <- as.factor(E(allcases.graph)$case)

V(allcases.graph)$label <- NA

E(allcases.graph)$color[E(allcases.graph)$case == 1] <- 'red'    
E(allcases.graph)$color[E(allcases.graph)$case == 2] <- 'blue'  
E(allcases.graph)$color[E(allcases.graph)$case == 4] <- 'red'  

# V(allcases.graph)$label.cex = 1.5

# plot(allcases.graph, layout = layout_with_fr,vertex.label.degree=0)
# plot(allcases.graph, layout = layout_with_fr, vertex.label.dist=2)
plot(allcases.graph, layout = layout_as_tree, vertex.label.dist=2, main="Case1 - TRUE and NA")
plot(allcases.graph, layout=layout_with_fr, vertex.label.dist=2, main="Case1 - TRUE and NA")
# plot(allcases.graph, layout = layout_with_fr, vertex.label.dist=1)








# start_time <- Sys.time()
# end_time <- Sys.time()
# end_time - start_time


