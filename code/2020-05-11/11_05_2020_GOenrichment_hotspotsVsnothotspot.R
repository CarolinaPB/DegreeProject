goframeData <- unique(genes_GO.bio[,.(GO.identifier, evidence, gene)])
gs <- getgeneset(goframeData)

# genes in hotspot
genesA <- genesinhotspot$gene
# genes not in hotspot
genesB <- names(phenotype[,2:ncol(phenotype)])[!names(phenotype[,2:ncol(phenotype)]) %in% genesinhotspot$gene]

# universe is genes involved in causality
universe <- names(phenotype[,2:ncol(phenotype)])

# get enrichment
hgCondA <- getenrichment(gs, universe = universe, interestinggenes = genesA, cond = T)
hgCondA.dt <- data.table(summary(hgCondA))

hgCondB <- getenrichment(gs, universe = universe, interestinggenes = genesB, cond = T)
hgCondB.dt <- data.table(summary(hgCondB))

# causal genes
tocombine.A <- data.table(-log10(hgCondA.dt$Pvalue), hgCondA.dt$Term)
tocombine.B <- data.table(-log10(hgCondB.dt$Pvalue), hgCondB.dt$Term)
combined <- merge(tocombine.A, tocombine.B, by="V2", all=T)
colnames(combined) <- c("term", "in hotspot", "not in hotspot")

combined[is.na(`in hotspot`)]$`in hotspot` <- 0
combined[is.na(`not in hotspot`)]$`not in hotspot` <- 0
plot_enrichment_heatmap(combined)
