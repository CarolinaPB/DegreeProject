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
goframeData <- unique(genes_GO.bioprocess[,.(Gene.ontologyAnnotations.ontologyTerm.identifier, Gene.ontologyAnnotations.evidence.code.code, Gene.secondaryIdentifier)])

# creating the sets of genes to use 
# the "interesting genes" - either genesA (the causal ones) or genesB (on the receiving end)
# the "universe" - in this case, only the genes that are involved in the causality
genesA <- unlist(unique(find.effects_TF[,geneA]))
genesB <- unlist(unique(find.effects_TF[,geneB]))

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


# universe is all the genes
universe.all <- names(phenotype[,2:ncol(phenotype)])

res.geneA.uniall <- getenrichment(gs, universe = universe.all, interestinggenes = genesA)
res.geneA.uniall.dt <- data.table(summary(res.geneA.uniall))
#htmlReport(res.geneA.parents, file=paste0(respath, "2020-02-13/GOenrichmentreport_genesA.html"))

# get enrichment for genesB
res.geneB.uniall <- getenrichment(gs, universe = universe.all, interestinggenes = genesB)
res.geneB.uniall.dt <- data.table(summary(res.geneB.uniall))


# check if I got the same results using the universe = all genes and universe = genes involved in causality
all.equal(res.geneA, res.geneA.uniall)
all.equal(res.geneB, res.geneB.uniall)


# get list with different graphs that represent relations between GO terms
termgrA <- termGraphs(res.geneA, use.terms = T, pvalue = 0.05)

#save all the graphs to a pdf
# pdf(file = "results/2020-02-19/termgraph_A_bioproc.pdf", onefile = T)
# for (i in 1:length(termgrA)){
#   plotGOTermGraph(termgrA[[i]], r = res.geneA, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=30)
# }
# dev.off()


# get list with different graphs that represent relations between GO terms
termgrB <- termGraphs(res.geneB, use.terms = T, pvalue = 0.05)
# pdf(file = "results/2020-02-19/termgraph_B_bioproc.pdf", onefile = T)
# for (i in 1:length(termgrB)){
#   plotGOTermGraph(termgrB[[i]], r = res.geneB, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=30)
# }
# dev.off()


#### CONDITIONAL HYPERGEO TEST ######

# for genesA as the genes of interest
paramsCondA <- GSEAGOHyperGParams(name="first try",
                                 geneSetCollection=gs,
                                 geneIds = genesA,
                                 universeGeneIds = universe,
                                 ontology = "BP",
                                 pvalueCutoff = 0.05,
                                 conditional = T,
                                 testDirection = "over")
hgCondA = hyperGTest(paramsCondA)
hgCondA.dt <- data.table(summary(hgCondA))


# for genesB as the genes of interest
paramsCondB <- GSEAGOHyperGParams(name="first try",
                                 geneSetCollection=gs,
                                 geneIds = genesB,
                                 universeGeneIds = universe,
                                 ontology = "BP",
                                 pvalueCutoff = 0.05,
                                 conditional = T,
                                 testDirection = "over")
hgCondB = hyperGTest(paramsCondB)
hgCondB.dt <- data.table(summary(hgCondB))

# to only keep GO terms with at least 100 gene annotations
summary(hgCondB, categorySize=100)



# GO terms that are marked significant by the standard hypergeo test, but not by the conditional test
stdIdsA = sigCategories(res.geneA)
condIdsA = sigCategories(hgCondA)
# num of GO terms that were not significant with the conditional hypergeo test
length(setdiff(stdIdsA, condIdsA))

stdIdsB = sigCategories(res.geneB)
condIdsB = sigCategories(hgCondB)
# num of GO terms that were not significant with the conditional hypergeo test
length(setdiff(stdIdsB, condIdsB)) 


# create HTML reports (tables with enrichment)
htmlReport(hgCondA, file="results/2020-02-19/hgCondA_htmlreport.md")
htmlReport(hgCondB, file="results/2020-02-19/hgCondB_htmlreport.md")

# for the causal genes:
# num of enriched GO terms
nrow(res.geneA.dt)

# categories not enriched after conditional hypergeo
res.geneA.dt[GOBPID %in% setdiff(stdIds, condIds)]



# get list with different graphs that represent relations between GO terms
termgrA.cond <- termGraphs(hgCondA, use.terms = T, pvalue = 0.05)
# save all the graphs to a pdf
pdf(file = "results/2020-02-19/termgraph_A_bioproc_cond.pdf", onefile = T)
# for (i in 1:length(termgrA.cond)){
#   plotGOTermGraph(termgrA.cond[[i]], r = hgCondA, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=30)
# }
# dev.off()


# get list with different graphs that represent relations between GO terms
termgrB.cond <- termGraphs(hgCondB, use.terms = T, pvalue = 0.05)
# pdf(file = "results/2020-02-19/termgraph_B_bioproc_cond.pdf", onefile = T)
# for (i in 1:length(termgrB.cond)){
#   plotGOTermGraph(termgrB.cond[[i]], r = hgCondB, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=30)
# }
# dev.off()


# check if the GO terms are still significant after the conditional hypergeo test
terms <- nodes(termgrA[[2]])
hgCondA.dt.sub <- data.table(summary(hgCondA, pval=0.5))#hgCondA.dt[,.(GOBPID,Term, Pvalue)]

hypergeo_compare <- merge(res.geneA.dt[, .(GOBPID, Term, Pvalue)], hgCondA.dt.sub[order(GOBPID)], by=c("GOBPID", "Term"), all.y=T)
setnames(hypergeo_compare, old=c("Pvalue.x", "Pvalue.y"), new=c("pval", "cond.pval"))
# Adds new column - if true, it means that that GO term's p-value changed
hypergeo_compare[, changed:=pval != cond.pval]

# these show that a term was removed with the conditional hypergeo
hypergeo_compare[GOBPID %in% terms]
res.geneA.dt[GOBPID %in% terms]




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


