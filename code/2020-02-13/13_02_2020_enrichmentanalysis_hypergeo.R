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

## table with GO terms + parent terms for each gene
ympy.evidence.parents <- fread(paste0(path, "results/2020-02-13/yeastmine_evidencecode_parents.gz"))
ympy.evidence.parents <- unique(ympy.evidence.parents)

goframeData <- unique(ympy.evidence.parents[,.(Gene.goAnnotation.ontologyTerm.identifier, Gene.goAnnotation.evidence.code.code, Gene.secondaryIdentifier)])


genesA <- unlist(unique(find.effects_TF[,geneA]))
genesB <- unlist(unique(find.effects_TF[,geneB]))
universe <- colnames(phenotype[,2:ncol(phenotype)])

gs <- getgeneset(goframeData)
res.geneA <- getenrichment(gs, universe = universe, interestinggenes = genesA)
res.geneA.dt <- data.table(summary(res.geneA))
#htmlReport(res.geneA.parents, file=paste0(respath, "2020-02-13/GOenrichmentreport_genesA.html"))

res.geneB <- getenrichment(gs, universe = universe, interestinggenes = genesB)
res.geneB.dt <- data.table(summary(res.geneB))


termgrA <- termGraphs(res.geneA, use.terms = T, pvalue = 0.05)
pdf(file = "results/2020-02-13/termgraph_A.pdf", onefile = T)
for (i in 1:length(termgrA)){
  plotGOTermGraph(termgrA[[i]], r = res.geneA, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=30)
}
dev.off()

termgrB <- termGraphs(res.geneB, use.terms = T, pvalue = 0.05)
pdf(file = "results/2020-02-13/termgraph_B.pdf", onefile = T)
for (i in 1:length(termgrB)){
  plotGOTermGraph(termgrB[[i]], r = res.geneB, add.counts = T, node.colors=c(sig="green", not="white"), max.nchar=30)
}
dev.off()


#### Conditional hypergeo test
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

res.geneA.dt[GOBPID %in% terms, .(GOBPID, Term, Pvalue)][order(GOBPID)] # standard
hgCond.dt.sub[GOBPID %in% terms][order(GOBPID)] # after conditional 
