##### LOAD FILES #####
library(data.table)
library(tidyr)
library(parallel)
library(igraph)
library(Hmisc)
library(corrplot)
library(dplyr)

# path <- "/home/carolpb/DegreeProject/" # use with UPPMAX
path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" # use with own computer
# path <- "/home/carolina/DegreeProject/" # use with diprotodon

phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))
#eqtl_results <- fread(paste0(path,"results/2019-12-11/eqtl_results_effect_eqtl-gene.gz"))
eqtl_results <- fread(paste0(path,"data/SI_Data_04_eQTL.csv"))
eqtl_results[cis=="VERDADEIRO"]$cis <- T
eqtl_results[cis=="FALSO"]$cis <- F


source("code/myfunctions.R")

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

genesB <- colnames(phenotype[,2:ncol(phenotype)])
effectsA_B.sepA_B <- create_ini_table(eqtl_results.sub, genesB, var.exp.lim)


# Final table with geneA and geneB and their corresponding eqtls (geneB might not have an eqtl)
effectsA_B.sepA_B

fwrite(effectsA_B.sepA_B,"results/2020-01-01/infoA_B.gz")


effectsA_B.sepA_B <- fread(paste0(path,"results/2020-01-07/infoA_B.gz"))


#######
### Find effect of eqtls from geneA in expression of geneB ####
eqtls.A <- unique(effectsA_B.sepA_B$eqtl.A)
genes.B <- unique(effectsA_B.sepA_B$geneB)
res.tot.eqtlA_geneB <- data.table(expand.grid(gene=genes.B, eqtl=eqtls.A))#, anv.res=NA))
res.eqtlA_geneB <- res.tot.eqtlA_geneB

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
res.eqtlB_geneA <- res.tot.eqtlB_geneA

# run the anova function in parallel -- effect of eqtlB on geneA
start_time <- Sys.time()
cl = makeCluster(detectCores() - 1, type="FORK")
res.eqtlB_geneA$anv.res <- parApply(cl=cl,res.eqtlB_geneA,1,effect_eqtl_gene, phenotype, genotype)
stopCluster(cl)
end_time <- Sys.time()
end_time - start_time


## merge anova results with the information table to create an effects table
# read in the anova results files
res.eqtlA_geneB <- fread(paste0(path, "results/2020-01-09/anova_eqtlA_geneB.gz"))
res.eqtlB_geneA <- fread(paste0(path, "results/2020-01-09/anova_eqtlB_geneA.gz"))

# results from effect of eqtlA on geneB
effects_table.eqtlA_geneB <- merge_after_anova(res.eqtlA_geneB, "B", "A", effectsA_B.sepA_B)

# results from effect of eqtlB on geneA
effects_table.anova <- merge_after_anova(res.eqtlB_geneA, gene.AB="A", eqtl.AB="B", effects_table.eqtlA_geneB )
setcolorder(effects_table.anova, c("geneA", "geneB", "eqtl.A", "eqtl.B", "cis.A", "cis.B"))

fwrite(effects_table.anova, paste0(path, "results/2020-01-09/09_01_2020_anovatable.gz"), na = NA)
effects_table.anova <- fread(paste0(path, "results/2020-01-09/09_01_2020_anovatable.gz"))





### correlation ####
cor_traits <- rcorr(as.matrix(phenotype[,2:ncol(phenotype)])) # to remove the sample names
save(cor_traits,file= "/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-01-08/full_cor.Rdata")
load(file = paste0(path, "results/2020-01-08/full_cor.Rdata"))

cor_traits.cor <- cor_traits$r
cor_traits.p <- cor_traits$P

my_cor_matr_flat <- flat_cor_mat(cor_traits.cor,cor_traits.p)
my_cor_matr_flat <- data.table(my_cor_matr_flat)

cor_matr <- my_cor_matr_flat[!duplicated(t(apply(my_cor_matr_flat, 1, sort))), ]

fwrite(cor_matr, paste0(path, "results/2020-01-10/cor_matr_unique.gz"))

cor_matr <- fread(paste0(path, "results/2020-01-10/cor_matr_unique.gz"))


# merge the correlation values with the table that has the anova results
setnames(cor_matr, old=c("row","column", "pval"), new=c("geneA","geneB", "cor.pval"))

effects_table.cor <- merge(effects_table.anova, cor_matr, by=c("geneA","geneB"), all.x=T)

fwrite(effects_table.cor, paste0(path, "results/2020-01-10/effectstable.gz"), na=NA)
effects_table.cor <- fread(paste0(path, "results/2020-01-10/effectstable.gz"))

######

# what I can change to get diff results:
# var.exp.lim
# snp.pval
# snp.pval.nsign
# corr.pval

sign_p <- c(1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
non_sign_p <- c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
params0 <- expand.grid(sign_p=sign_p, non_sign_p=non_sign_p)







find.effects <- effects_table.cor[cor.pval < corr.pval & cis.A ==T & cis.B==T & geneA!=geneB]

find.effects <- find.effects_fun(find.effects, snp.pval, snp.pval.nsign)

####
# find.effects[(find.effects$`A->B`==T & find.effects$`B->A`==F) | (find.effects$`A->B`==F & find.effects$`B->A`==T)]


find.effects_TF.1 <- find.effects[find.effects$`A->B`==T & find.effects$`B->A`==F, .(geneA, geneB, eqtl.A, eqtl.B, `A->B`, `B->A`)]
find.effects_TF.2 <- rbind(find.effects_TF.1, find.effects[find.effects$`A->B`==F & find.effects$`B->A`==T, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
find.effects_TF <- unique(find.effects_TF.2)


# when one direction is true and the other is false
find.effects_TF.numpairs <- find.effects_TF[(`A->B`==T & `B->A`==F), .N, by=.(geneA, geneB)]

numpairs.table <- data.table(N = unique(find.effects_TF.numpairs$N), numpairs =NA)
for (i in 1:nrow(numpairs.table)){
  numpairs.table$numpairs[i] <- nrow(find.effects_TF.numpairs[N==numpairs.table$N[i]])
}

# plot the number of times a gene pair appears
bp <- barplot(numpairs.table$numpairs, numpairs.table$N, names.arg=unique(find.effects_TF.numpairs$N), 
              width = 0.5, space=0.2, legend.text = F, ylim = c(0,max(numpairs.table$numpairs)+2000),
              main = "Number of times gene pairs appear \n (A->B = T and B->A = F)", xlab = "# times a gene pair appears", 
              ylab = "# of gene pairs")
text(bp,numpairs.table$numpairs, labels=numpairs.table$numpairs, cex=1, pos=3)

# most gene pairs appear 1 time

find.effects_TF.geneeqtl.A <- find.effects_TF[(`A->B`==T & `B->A`==F), .N, by=.(geneA, eqtl.A)]
find.effects_TF.geneeqtl.B <- find.effects_TF[(`A->B`==T & `B->A`==F), .N, by=.(geneB, eqtl.B)]

find.effects_TF.geneeqtl.A.plot <- find.effects_TF.geneeqtl.A %>% unite(gene_eqtlA, geneA, eqtl.A, sep = "__")
#find.effects_TF.geneeqtl.B.plot <- find.effects_TF.geneeqtl.B %>% unite(gene_eqtlB, geneB, eqtl.B, sep = "__")

plot(as.factor(find.effects_TF.geneeqtl.A.plot$gene_eqtlA), find.effects_TF.geneeqtl.A.plot$N, 
     xlab="gene-eqtl pairs", ylab="# times each pair appears", axes=FALSE, 
     main="Num times of times each gene-eqtl pair appears \n (A->B = T and B->A = F)")
Axis(side=2, labels=T)


# freq of gene-eqtl pairs
numpairs.table2 <- data.table(N = unique(find.effects_TF.geneeqtl.A.plot$N), numpairs =NA)
for (i in 1:nrow(numpairs.table2)){
  numpairs.table2$numpairs[i] <- nrow(find.effects_TF.geneeqtl.A[N==numpairs.table2$N[i]])
}
numpairs.table2 <- numpairs.table2[order(N)]

# all values individually
bp2 <- barplot(numpairs.table2$numpairs, numpairs.table2$N, names.arg=unique(numpairs.table2$N), 
              width = 0.5, space=0.2, legend.text = F, ylim = c(0,max(numpairs.table2$numpairs)+100),
              main = "Number of times gene-eqtl pairs appear \n (A->B = T and B->A = F)", xlab = "# times a gene-eqtl pair appears", 
              ylab = "Freq of gene-eqtl pairs", cex.names=0.8)

# by binning values 
hist(find.effects_TF.geneeqtl.A.plot$N, main = "Frequency of gene-eqtl pairs", xlab = "# times gene-eqtl pairs appear",labels = T)
# the majority of gene-eqtl pairs appear between 0-50 times

















#####

effects <- find.effects

# Define cases
# one gene affects the other but we don't know the action of the other
case1 <- effects[(is.na(effects$`A->B`) & effects$`B->A` == TRUE) | (is.na(effects$`B->A`) & effects$`A->B` == TRUE)]
# one gene affects the other but is not affected by the other gene
case2 <- effects[(effects$`A->B` == TRUE & effects$`B->A` == FALSE) | (effects$`A->B` == FALSE & effects$`B->A` == TRUE)]
# we don't know if the genes affect each other
case3 <- effects[is.na(effects$`A->B`) & is.na(effects$`B->A`)]
# a gene affects the other but it's also affected by that other gene
case4 <- effects[(effects$`A->B` == TRUE & effects$`B->A` == TRUE)]
# no gene affects the other
case5 <- effects[(effects$`A->B` == FALSE & effects$`B->A` == FALSE)]




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




#### Case2 plot
allcases.graph <- graph_from_data_frame(allcases.links[allcases.links$case == 2],directed=TRUE)

E(allcases.graph)$color <- as.factor(E(allcases.graph)$case)

V(allcases.graph)$label <- NA

E(allcases.graph)$color[E(allcases.graph)$case == 1] <- 'red'    
E(allcases.graph)$color[E(allcases.graph)$case == 2] <- 'blue'  
E(allcases.graph)$color[E(allcases.graph)$case == 4] <- 'red'  

# V(allcases.graph)$label.cex = 1.5

# plot(allcases.graph, layout = layout_with_fr,vertex.label.degree=0)
# plot(allcases.graph, layout = layout_with_fr, vertex.label.dist=2)
plot(allcases.graph, layout = layout_as_tree, vertex.label.dist=2, main="Case2 - TRUE and FALSE")
plot(allcases.graph, layout=layout_with_graphopt, vertex.label.dist=2, main="Case2 - TRUE and FALSE")
########
find.effects

# find.effects[,.(count= .N, count=sum(eqtl.A), by=geneA)]

### TODO table with summary numbers
### TODO plot with diff params
# ex: n gene pairs with A->B TRUE and B->A FALSE
# other ex:
# The summary table shows
# * number of unique traits
# * num of unique snps
# * num of A->B
# * num of B->A
# * num of A->B that we can’t say if it’s true or false
# * num of B->A that we can’t say if it’s true or false
# * num of unique gene pairs with A->B
# * num of unique gene pairs with B->A
# * num of unique gene pairs with A not affecting B
# * num of unique gene pairs with B not affecting A
# * num of unique gene pairs with undertermined A->B
# * num of unique gene pairs with undertermined B->A

res <- find.effects
tab <- data.table(matrix(ncol=16, nrow=1))
names(tab) <- c("geneA", "geneB", "eqtl.A", "eqtl.B","A->B", "B->A", "maybeAB","maybeBA", "A->B.gpairs", "B->A.gpairs", "~A->B.gpairs","~B->A.gpairs", "A->Bmaybe.gpairs","B->Amaybe.gpairs", "B->A.gpairs_TF", "A->B.gpairs_TF")
tab$geneA <- length(unique(res$geneA)) # num of unique genes A
tab$geneB <- length(unique(res$geneB)) # num of unique genes B
tab$eqtl.A <- length(unique(res$eqtl.A)) # num of unique eqtlA
tab$eqtl.B <- length(unique(res$eqtl.B)) # num of unique eqtlB
tab$maybeBA <- nrow(res[is.na(res$`B->A`)]) # num of gene-eqtl pairs for which we can't say
tab$maybeAB <- nrow(res[is.na(res$`A->B`)]) # num of gene-eqtl pairs for which we can't say
tab$`A->B` <- nrow(res[res$`A->B`==TRUE])
tab$`B->A` <- nrow(res[res$`B->A`==TRUE])
tab$`A->B.gpairs` <- nrow(res[res$`A->B`==TRUE, .N, by=.(geneA, geneB)])  # number of gene pairs for which A->B is true
tab$`B->A.gpairs` <- nrow(res[res$`B->A`==TRUE, .N, by=.(geneA, geneB)])  # number of gene pairs for which B->A is true
tab$`A->Bmaybe.gpairs`<- nrow(res[is.na(res$`A->B`), .N, by=.(geneA, geneB)]) # number of gene pairs for which A->B is maybe
tab$`B->Amaybe.gpairs`<- nrow(res[is.na(res$`B->A`), .N, by=.(geneA, geneB)]) # number of gene pairs for which B->A is maybe
tab$`~A->B.gpairs` <- nrow(res[res$`A->B`==FALSE, .N, by=.(geneA, geneB)])  # number of gene pairs for which A->B is false
tab$`~B->A.gpairs` <- nrow(res[res$`B->A`==FALSE, .N, by=.(geneA, geneB)])  # number of gene pairs for which B->A is false
tab$`A->B.gpairs_TF` <- nrow(res[res$`B->A`==FALSE & res$`A->B`==TRUE, .N, by=.(geneA, geneB)])  # number of gene pairs for which B->A is false
tab$`B->A.gpairs_TF` <- nrow(res[res$`B->A`==TRUE & res$`A->B`==FALSE, .N, by=.(geneA, geneB)])  # number of gene pairs for which B->A is false

t(tab)

A_to_B <- res[res$`B->A`==FALSE & res$`A->B`==TRUE, .N, by=.(geneA, geneB)]
B_to_A <- res[res$`B->A`==TRUE & res$`A->B`==FALSE, .N, by=.(geneA, geneB)]




# Find number of each geno for each SNP ######
snp_ids <- colnames(genotype)[2:length(colnames(genotype))]
# nSNPs <- length(snp_ids)

ngenos.dt <- data.table(snp_ids)

system.time({
  cl = makeCluster(detectCores() - 1, type="FORK")
  ngenos.dt$nums <- parApply(cl=cl,ngenos.dt,1,get_num_genos, genotype)
  stopCluster(cl)
})

numgenos <- ngenos.dt %>% separate(nums, c("n1", "n-1"), sep="__")

fwrite(numgenos, paste0(path, "results/2020-01-08/08_01_2020_numgenos.gz"))

numgenos <- fread(paste0(path, "results/2020-01-08/08_01_2020_numgenos.gz"))
# all genos are represented by a large amount of samples ####