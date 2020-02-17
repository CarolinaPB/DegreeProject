##### LOAD FILES #####
library(data.table)
library(tidyr)
library(parallel)
library(igraph)
library(dplyr)

# path <- "/home/carolpb/DegreeProject/" # use with UPPMAX
path <- "/Users/Carolina/Documents/GitHub/DegreeProject/"
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"# use with own computer
# path <- "/home/carolina/DegreeProject/" # use with diprotodon

source("code/myfunctions.R")

##### PARAMETERS #####
var.exp.lim <- 0.1

# nSNPs <- length(colnames(genotype))-1
# nGenes <- length(colnames(phenotype))-1

nSNPs <- 42052
nGenes <- 5720

snp.pval <- 0.01
snp.pval.nsign <- as.numeric(1e-5)

corr.pval <- 0.05/choose(nGenes,2)


######
effects_table.cor <- fread(paste0(path, "results/2020-01-10/effectstable.gz"))

find.effects <- effects_table.cor[cor.pval < corr.pval & cis.A ==T & cis.B==T & geneA!=geneB]

find.effects <- find.effects_fun(find.effects, snp.pval, snp.pval.nsign)

####
find.effects_TF.1 <- find.effects[find.effects$`A->B`==T & find.effects$`B->A`==F, .(geneA, geneB, eqtl.A, eqtl.B, `A->B`, `B->A`)]
find.effects_TF.2 <- rbind(find.effects_TF.1, find.effects[find.effects$`A->B`==F & find.effects$`B->A`==T, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
find.effects_TF <- unique(find.effects_TF.2)

#Number of times each gene is the causal one or is on the receiving end
ngeneA <- find.effects_TF[find.effects_TF$`A->B`==T, .N, by=geneA]
ngeneB <- find.effects_TF[find.effects_TF$`B->A`==F, .N, by=geneB]
ngeneA_B <- merge(ngeneA, ngeneB, by.x = "geneA", by.y = "geneB", all=T)
colnames(ngeneA_B) <- c("gene", "A","B")
ngeneA_B[order(-A)]


# when one direction is true and the other is false
find.effects_TF.numpairs <- find.effects_TF[(`A->B`==T & `B->A`==F), .N, by=.(geneA, geneB)]

numpairs.table <- data.table(N = unique(find.effects_TF.numpairs$N), numpairs =NA)
for (i in 1:nrow(numpairs.table)){
  numpairs.table$numpairs[i] <- nrow(find.effects_TF.numpairs[N==numpairs.table$N[i]])
}

# plot the number of times a gene pair appears
bp <- barplot(numpairs.table[order(-numpairs)]$numpairs, numpairs.table[order(-numpairs)]$N, names.arg=unique(find.effects_TF.numpairs[order(N)]$N), 
              width = 0.5, space=0.2, legend.text = F, ylim = c(0,max(numpairs.table$numpairs)+2500),
              main = "Number of times gene pairs appear \n (A->B = T and B->A = F)", xlab = "# times a gene pair appears", 
              ylab = "# of gene pairs")
text(bp,numpairs.table[order(-numpairs)]$numpairs, labels=numpairs.table[order(-numpairs)]$numpairs, cex=1, pos=3)


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
               ylab = "# gene-eqtl pairs", cex.names=0.8)

# by binning values 
hist(find.effects_TF.geneeqtl.A.plot$N, main = "# gene-eqtl pairs", xlab = "# times gene-eqtl pairs appear",labels = T)
# the majority of gene-eqtl pairs appear between 0-50 times


##### network plotting ####
plot_TF <- graph_from_edgelist(as.matrix(find.effects_TF[,.(geneA, geneB)]),directed=TRUE)
# components(plot_TF)

memb <- components(plot_TF)$membership
nodes <- data.table(id = names(memb),
                    group_id = memb)
nodes <- nodes[order(nodes$id), ]
nodes

plot(plot_TF, layout=layout_with_gem,edge.arrow.size=.5, vertex.label=NA)






find.effects_TNA.1 <- find.effects[find.effects$`A->B`==T & is.na(find.effects$`B->A`), .(geneA, geneB, eqtl.A, eqtl.B, `A->B`, `B->A`)]
find.effects_TNA.2 <- rbind(find.effects_TNA.1, find.effects[is.na(find.effects$`A->B`) & find.effects$`B->A`==T, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
find.effects_TNA <- unique(find.effects_TNA.2)

plot_TNA <- graph_from_edgelist(as.matrix(find.effects_TNA[,.(geneA, geneB)]),directed=TRUE)
# components(plot_TF)

memb <- components(plot_TNA)$membership
nodes <- data.table(id = names(memb),
                    group_id = memb)
nodes <- nodes[order(nodes$id), ]
nodes

V(plot_TNA)$size <- 5
plot(plot_TNA, layout=layout_with_gem,edge.arrow.size=.5, vertex.label=NA, edge.curved=.1)

