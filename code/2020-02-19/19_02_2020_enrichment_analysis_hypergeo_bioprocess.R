library("GSEABase")
library("GOstats")
library(parallel)
library(data.table)
# 
path <- "/Users/Carolina/Documents/GitHub/DegreeProject/"
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"

# path <- "/home/carolpb/DegreeProject/" # uppmax
# respath <- "/proj/snic2019-8-367/private/carol/results/" # uppmax

source(paste0(path, "code/myfunctions.R"))

phenotype <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/data/SI_Data_01_expressionValues.txt")

find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))

## table with GO terms + parent terms for each gene + gene evidence code
genes_GO.table <- fread(paste0(path, "results/2020-02-19/genelistwithGOterm_allnamespaces.tsv"))
genes_GO.table <- unique(genes_GO.table)
genes_GO.bioprocess <- genes_GO.table[Gene.ontologyAnnotations.ontologyTerm.namespace=="biological_process"]
genes_GO.bioprocess[,Gene.ontologyAnnotations.ontologyTerm.namespace:=NULL]



# table with necessary columns for creating a geneset and calculating enrichment
# NOT USING PARENTS INFO - USING THE "NORMAL" INFO
goframeData <- unique(genes_GO.bioprocess[,.(Gene.ontologyAnnotations.ontologyTerm.identifier, Gene.ontologyAnnotations.evidence.code.code, Gene.secondaryIdentifier)])

# creating the sets of genes to use 
# the "interesting genes" - either genesA (the causal ones) or genesB (on the receiving end)
# the "universe" - in this case, all the genes in our dataset
genesA <- unlist(unique(find.effects_TF[,geneA]))
genesB <- unlist(unique(find.effects_TF[,geneB]))
# universe <- colnames(phenotype[,2:ncol(phenotype)])
universe <- unique(c(genesA, genesB))
# create geneset
gs <- getgeneset(goframeData)
# get enrichment for genesA
res.geneA <- getenrichment(gs, universe = universe, interestinggenes = genesA)
res.geneA.dt <- data.table(summary(res.geneA))
#htmlReport(res.geneA.parents, file=paste0(respath, "2020-02-13/GOenrichmentreport_genesA.html"))

# get enrichment for genesB
res.geneB <- getenrichment(gs, universe = universe, interestinggenes = genesB)
res.geneB.dt <- data.table(summary(res.geneB))

# get list with different graphs that represent relations between GO terms
termgrA <- termGraphs(res.geneA, use.terms = T, pvalue = 0.05)

#save all the graphs to a pdf
pdf(file = "results/2020-02-19/termgraph_A_bioproc.pdf", onefile = T)
for (i in 1:length(termgrA)){
  plotGOTermGraph(termgrA[[i]], r = res.geneA, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=30)
}
dev.off()


# get list with different graphs that represent relations between GO terms
termgrB <- termGraphs(res.geneB, use.terms = T, pvalue = 0.05)
pdf(file = "results/2020-02-19/termgraph_B_bioproc.pdf", onefile = T)
for (i in 1:length(termgrB)){
  plotGOTermGraph(termgrB[[i]], r = res.geneB, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=30)
}
dev.off()


#### CONDITIONAL HYPERGEO TEST ######

# for genesA as the genes of interest
paramsCond <- GSEAGOHyperGParams(name="first try",
                                 geneSetCollection=gs,
                                 geneIds = genesA,
                                 universeGeneIds = universe,
                                 ontology = "BP",
                                 pvalueCutoff = 0.05,
                                 conditional = T,
                                 testDirection = "over")
hgCond = hyperGTest(paramsCond)
hgCond.dt <- data.table(summary(hgCond))


# GO terms that are marked significant by the standard hypergeo test, but not by the conditional test
stdIds = sigCategories(res.geneA)
condIds = sigCategories(hgCond)
setdiff(stdIds, condIds)


# check if the GO terms are still significant after the conditional hypergeo test
terms <- nodes(termgrA[[2]])
hgCond.dt.sub <- data.table(summary(hgCond, pvalue=0.5)[,c("GOBPID","Term", "Pvalue")])

hypergeo_compare <- merge(res.geneA.dt[GOBPID %in% terms, .(GOBPID, Term, Pvalue)], hgCond.dt.sub[GOBPID %in% terms][order(GOBPID)], by=c("GOBPID", "Term"), all=T)
setnames(hypergeo_compare, old=c("Pvalue.x", "Pvalue.y"), new=c("pval", "cond.pval"))
# Adds new column - if true, it means that that GO term's p-value changed
hypergeo_compare[, changed:=pval != cond.pval]

hypergeo_compare

# "cellular response to DNA damage stimulus" is no longer significant -  does it mean that we can only say that only the general term DNA repair is significant?

# check if the GO terms are still significant after the conditional hypergeo test
# first node group
terms.1 <- nodes(termgrA[[1]])

res.geneA.dt[GOBPID %in% terms.1, .(GOBPID, Term, Pvalue)][order(GOBPID)] # standard
hgCond.dt.sub[GOBPID %in% terms.1][order(GOBPID)] # after conditional 


hypergeo_compare.1 <- merge(res.geneA.dt[GOBPID %in% terms.1, .(GOBPID, Term, Pvalue)], hgCond.dt.sub[GOBPID %in% terms.1][order(GOBPID)], by=c("GOBPID", "Term"), all=T)
setnames(hypergeo_compare.1, old=c("Pvalue.x", "Pvalue.y"), new=c("pval", "cond.pval"))

hypergeo_compare.1[, changed:=pval != cond.pval]
hypergeo_compare.1

# terms that changed but are still significant at pval < 0.05
hypergeo_compare[pval<cond.pval & cond.pval < 0.05]




termgrA_cond <- termGraphs(hgCond, use.terms = T, pvalue = 0.05)
# pdf(file = "results/2020-02-13/termgraph_A_condhypergeo.pdf", onefile = T)
# for (i in 1:length(termgrA_cond)){
#   plotGOTermGraph(termgrA_cond[[i]], r = hgCond, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=200)
# }
# dev.off()

