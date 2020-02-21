eqtl_results <- fread("data/SI_Data_04_eQTL.csv")


find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))

genes_coord <- eqtl_results[,.(gene, genePlotCoordinate)]

plot(eqtl_results$markerPlotCoordinate, genes_coord$genePlotCoordinate, pch=".")

phenotype <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/data/SI_Data_01_expressionValues.txt")
fwrite(data.table(colnames(phenotype[,2:ncol(phenotype)])), "results/2020-02-20/allgenes_list.txt")

genepos <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-20/gene_pos.gz")
plot(genepos[,.(start, end)], pch=".")

causal.pos.A <- merge(find.effects_TF, genepos, by.x="geneA", by.y="gene", all.x=T)
causal.pos.B <- merge(causal.pos.A, genepos, by.x="geneB", by.y="gene", all.x=T)
plot(causal.pos.B$start.x, causal.pos.B$start.y, pch=".")

numlinks_from_geneA <- unique(find.effects_TF[,.(geneA, geneB)])[, .(l_out=.N), by=geneA]
numlinks_from_geneB <- unique(find.effects_TF[,.(geneA, geneB)])[, .(l_in=.N), by=geneB]
# links_count <- merge(numlinks_from_geneA, numlinks_from_geneB, by="geneB")
# setcolorder(links_count, c("geneA", "geneB", "countA", "countB"))

genes <- data.table(genes=unique(c(find.effects_TF$geneA, find.effects_TF$geneB)))
merge_lin <- merge(genes, numlinks_from_geneA, by.="genes", by.y="geneA", all.x=T)
genes_nlinks <- merge(merge_lin, numlinks_from_geneB, by.x="genes", by.y="geneB", all.x=T)

causalgenes <- genes_nlinks[l_out > 10]
affectedgenes <- genes_nlinks[l_out <2 & !is.na(l_in)]


# causalgenes <- unique(links_count[countA>10, .(geneA, countA)])
# affectedgenes <- unique(links_count[, .(geneA, countA)])

plot(causal.pos.B[geneA %in% causalgenes$genes]$start.x, causal.pos.B[geneA %in% affectedgenes$genes]$start.y, pch=".")
plot(genepos[,.(start, end)], pch=".")
lines(genepos[gene %in% causalgenes$genes]$start, genepos[gene %in% affectedgenes$genes]$start)


genes_nlinks[l_out >40]
genes_nlinks[is.na(l_out)]
