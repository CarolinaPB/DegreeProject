setwd("/Users/Carolina/Documents/GitHub/DegreeProject/summarizing/")
library(data.table)
load("data/counts.RData")

# 1 if the value if the gene is absent in that individual (if the fiels value is zero), 0 if it's present
getzeros <- apply(counts$pheno, 1, function(x) {ifelse(x==0, 1, 0)})
# sum the number of individuals that don't have that gene
howmany <- colSums(getzeros)
# how many samples don't have a certain gene
count_missinggene <- data.table(gene=names(howmany), nsamples_missing =howmany)
count_missinggene <- count_missinggene[order(-nsamples_missing)]



find.effects_TF <- fread("results/findeffects_TF_newparams.gz")

test <- merge(unique(find.effects_TF[,.(geneA, eqtl.A)]), count_missinggene, by.x="geneA", by.y="gene")
test <- test[order(-nsamples_missing)]
