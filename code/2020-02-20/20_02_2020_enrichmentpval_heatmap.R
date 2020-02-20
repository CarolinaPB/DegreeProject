library(pheatmap)
library(data.table)
library("GSEABase")
library("GOstats")

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



tocombine.A <- data.table(-log10(hgCondA.dt$Pvalue), hgCondA.dt$Term)
tocombine.B <- data.table(-log(hgCondB.dt$Pvalue), hgCondB.dt$Term)
combined <- merge(tocombine.A, tocombine.B, by="V2", all=T)
colnames(combined) <- c("term", "genesA", "genesB")

Causal <- combined$genesA
Affected <- combined$genesB

test2 <- rbind(Causal, Affected)
colnames(test2) <- combined$term
test2[is.na(test2)] <- as.double(NA)
test2[is.na(test2)] <- 0

# pdf("results/results_figures/heatmap_enrichmentpvals.pdf")
ttest <- t(test2)




# pdf("results/results_figures/heatmap_enrichmentpvals.pdf")
paletteLength <- 50
myColor <- colorRampPalette(c("white", "yellow", "red"))(paletteLength)

myBreaks <- c(seq(min(ttest), max(ttest)/paletteLength-0.000001, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(ttest)/paletteLength, max(ttest), length.out=floor(paletteLength/2)))
pheatmap(ttest, display_numbers = matrix(ifelse(ttest == 0, "*******", ""), nrow(ttest)), 
         fontsize_row=3, cellwidth=70, cluster_col=F, cluster_rows=T, number_color="White", 
         fontsize = 5, treeheight_row = 0, color=myColor, breaks=myBreaks, border_color="white", 
         main="Heatmap enrichment p-value")
dev.off()
# 