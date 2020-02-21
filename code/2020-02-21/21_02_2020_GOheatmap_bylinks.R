library(pheatmap)
library(data.table)
library("GSEABase")
library("GOstats")

path <- "/Users/Carolina/Documents/GitHub/DegreeProject/"
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"

# path <- "/home/carolpb/DegreeProject/" # uppmax
# respath <- "/proj/snic2019-8-367/private/carol/results/" # uppmax

source(paste0(path, "code/myfunctions.R"))

phenotype <-fread("/Users/Carolina/Documents/GitHub/DegreeProject/data/SI_Data_01_expressionValues.txt")

find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))

## table with GO terms + parent terms for each gene + gene evidence code
genes_GO.table <-fread(paste0(path,"results/2020-02-19/genelistwithGOterm_allnamespaces.tsv"))
genes_GO.table <- unique(genes_GO.table)
genes_GO.bioprocess <- genes_GO.table[Gene.ontologyAnnotations.ontologyTerm.namespace == "biological_process"]
genes_GO.bioprocess[, Gene.ontologyAnnotations.ontologyTerm.namespace := NULL]



# table with necessary columns for creating a geneset and calculating enrichment
goframeData <- unique(genes_GO.bioprocess[, .(Gene.ontologyAnnotations.ontologyTerm.identifier,Gene.ontologyAnnotations.evidence.code.code,Gene.secondaryIdentifier)])


numlinks_from_geneA <- unique(find.effects_TF[, .(geneA, geneB)])[, .(l_out = .N), by = geneA]
numlinks_from_geneB <- unique(find.effects_TF[, .(geneA, geneB)])[, .(l_in = .N), by = geneB]
# links_count <- merge(numlinks_from_geneA, numlinks_from_geneB, by="geneB")
# setcolorder(links_count, c("geneA", "geneB", "countA", "countB"))

# table with how many links go in and out of a gene
genes <-data.table(genes = unique(c(find.effects_TF$geneA, find.effects_TF$geneB)))
merge_lin <-
  merge( genes, numlinks_from_geneA, by. = "genes", by.y = "geneA", all.x = T)
genes_nlinks <-merge( merge_lin, numlinks_from_geneB, by.x = "genes", by.y = "geneB", all.x = T)

# scatterplot of link distribution
# plot(genes_nlinks[, .(l_out, l_in)], main = "How many link each gene has going in or out", pch =".")
# plot(genes_nlinks[, .(l_out, l_in)], xlim = c(0, 200), pch = ".")

causalgenes <- genes_nlinks[l_out > 40]
affectedgenes <- genes_nlinks[is.na(l_out) & !is.na(l_in)]




# creating the sets of genes to use
# the "interesting genes" - either genesA (the causal ones) or genesB (on the receiving end)
# the "universe" - in this case, only the genes that are involved in the causality
genesA <- causalgenes$genes
genesB <- affectedgenes$genes

# universe <- unique(c(genesA, genesB))
universe <- unlist(unique(c(find.effects_TF[, geneA], find.effects_TF[, geneB])))
# universe2 <- names(phenotype[,2:ncol(phenotype)])

# create geneset
gs <- getgeneset(goframeData)
# get enrichment for genesA

hgCondA <-
  getenrichment(
    geneset = gs,
    universe = universe,
    interestinggenes = genesA,
    cond = TRUE
  )
hgCondA.dt <- data.table(summary(hgCondA))


# for genesB as the genes of interest
hgCondB = getenrichment(
  geneset = gs,
  universe = universe,
  interestinggenes = genesB,
  cond = TRUE
)
hgCondB.dt <- data.table(summary(hgCondB))


# table with all enriched GO terms found for both the causal and the affected genes with the corresponding enrichment p-val
tocombine.A <-
  data.table(-log10(hgCondA.dt$Pvalue), hgCondA.dt$Term)
tocombine.B <-
  data.table(-log10(hgCondB.dt$Pvalue), hgCondB.dt$Term)
combined <- merge(tocombine.A, tocombine.B, by = "V2", all = T)
colnames(combined) <- c("term", "genesA", "genesB")
combined[is.na(genesA)]$genesA <- 0
combined[is.na(genesB)]$genesB <- 0


# plot GO enrichment heatmap
plot_enrichment_heatmap(combined)
