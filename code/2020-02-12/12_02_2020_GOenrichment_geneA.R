######################################
# Getting enrichment score for geneA #
######################################

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

find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))

# python 10_02_2020_yeastminepy_go_evidencecode.py > /Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-10/yeastmine/yeastmine_evidencecode.txt
ympy.evidence <- fread(paste0(respath, "2020-02-10/yeastmine/yeastmine_evidencecode.txt"))


goframeData <- unique(ympy.evidence[,.(GO, evidencecode, gene)])

gs <- getgeneset(goframeData)
res.geneA <- getenrichment(gs, "geneA")
res.geneA.dt <- data.table(summary(res.geneA))

res.geneB <- getenrichment(gs, "geneB")
res.geneB.dt <- data.table(summary(res.geneB))

# which Terms are shared between genesA and genesB
res.geneB.dt[Term %in% res.geneA.dt$Term]
# to access the functions to extract info from the result of getenrichment/hyperGTest ?geneCounts




# TODO
# merge findeffects with GO for both geneA and B
ympy.GO.all <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastminepy/yeastmine_ontology_py.txt")
ympy.GO <- unique(ympy.GO.all[,.(gene, ontologyTerm.identifier, ontologyTerm.name)])
ympy.GO.parents <- unique(ympy.GO.all[,.(gene, ontologyTerm.parents.identifier, ontologyTerm.parents.name)])

find.effects.GO.geneA <- merge(find.effects_TF, ympy.GO, by.x="geneA", by.y="gene",all.x=T, allow.cartesian=T)
find.effects.GO <- merge(find.effects.GO.geneA, ympy.GO, by.x="geneB", by.y="gene",all.x=T, allow.cartesian=T)
find.effects.GO.sub <- unique(find.effects.GO[,.(geneA, geneB, ontologyTerm.identifier.x, ontologyTerm.identifier.y, ontologyTerm.name.x, ontologyTerm.name.y)])
setnames(find.effects.GO.sub, old=c("ontologyTerm.identifier.x", "ontologyTerm.identifier.y", "ontologyTerm.name.x", "ontologyTerm.name.y"), new=c("GO.A", "GO.B", "Term.A", "Term.B"))



unique(find.effects.GO.sub[, .(geneA, GO.A, Term.A)])

genes <- unlist(unique(find.effects_TF[,geneA]))
termGraphs(res.geneA, id = geneIds(res.geneA), use.terms = T)
inducedTermGraph(res.geneA, geneIds(res.geneA), children = T)

t <- termGraphs(res.geneA, use.terms = T, pvalue = 0.05)
pdf(file = "results/2020-02-13/termgraph1.pdf")
plotGOTermGraph(t$`1`, r = res.geneA, add.counts = T, node.colors=c(sig="green", not="white"))
dev.off()


# Use parent GO term instead of the specific GO code
find.effects.GOparents.geneA <- merge(find.effects_TF, ympy.GO.parents, by.x="geneA", by.y="gene",all.x=T, allow.cartesian=T)
find.effects.GOparents <- merge(find.effects.GOparents.geneA, ympy.GO.parents, by.x="geneB", by.y="gene",all.x=T, allow.cartesian=T)
# find.effects.GOparents.sub <- unique(find.effects.GOparents[,.(geneA, geneB, ontologyTerm.parents.identifier.x, ontologyTerm.parents.identifier.y, ontologyTerm.parents.name.x, ontologyTerm.parents.name.y)])
# setnames(find.effects.GO.sub, old=c("ontologyTerm.identifier.x", "ontologyTerm.identifier.y", "ontologyTerm.name.x", "ontologyTerm.name.y"), new=c("GO.A", "GO.B", "Term.A", "Term.B"))

ympy.evidence.parents <- fread(paste0(path, "results/2020-02-13/yeastmine_evidencecode_parents.gz"))
ympy.evidence.parents <- unique(ympy.evidence.parents)
goframeData <- unique(ympy.evidence.parents[,.(Gene.goAnnotation.ontologyTerm.parents.identifier, Gene.goAnnotation.evidence.code.code, Gene.secondaryIdentifier)])

phenotype <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/data/SI_Data_01_expressionValues.txt")


genes <- unlist(unique(find.effects_TF[,geneA]))
universe <- colnames(phenotype[,2:ncol(phenotype)])

gs.parents <- getgeneset(goframeData)
res.geneA.parents <- getenrichment(gs.parents, universe = universe, interestinggenes = genes)
res.geneA.dt.parents <- data.table(summary(res.geneA.parents))

res.geneB.parents <- getenrichment(gs.parents, "geneB")
res.geneB.dt.parents <- data.table(summary(res.geneB.parents))


t.parents <- termGraphs(res.geneA, use.terms = T, pvalue = 0.05)
pdf(file = "results/2020-02-13/termgraph1_parents.pdf")
plotGOTermGraph(t$`1`, r = res.geneA, add.counts = T, node.colors=c(sig="green", not="white"))
dev.off()


# which Terms are shared between genesA and genesB
res.geneB.dt[Term %in% res.geneA.dt$Term]
# to access the functions to extract info from the result of getenrichment/hyperGTest ?geneCounts


#TODO
# for every row, keep rows that have term.B = term and term.A = parent term

test <- find.effects.GO.sub[GO.A %in% ympy.GO.all[ontologyTerm.identifier %in% GO.B]$ontologyTerm.parents.identifier]
dim(find.effects.GO.sub)
dim(test)


find.effects.GO.sub.geneA <- unique(find.effects.GO.sub[, .(geneA,geneB, GO.A, Term.A)])
find.effects.GO.sub.geneA[,("sum"):=.N, by=c(Term.A)]
