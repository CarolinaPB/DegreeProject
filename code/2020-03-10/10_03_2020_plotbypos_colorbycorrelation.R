library(data.table)

genepos <- fread("results/gene_pos.gz")
if (!exists("find.effects_TF")){
  find.effects_TF <- fread("results/findeffects_TF_newparams.gz")
}
find.effects <- fread("results/findeffects_all_newparams.gz")

# Add chomosome and start position to each gene
causal.pos.A <- merge(find.effects_TF, genepos, by.x="geneA", by.y="gene", all.x=T)
colnames(causal.pos.A) <- c("geneA", "geneB", "eqtl.A", "eqtl.B", "A->B", "B->A", "strand.A", "start.A","end.A", "chr.A")
causal.pos.B <- merge(causal.pos.A, genepos, by.x="geneB", by.y="gene", all.x=T)


# Keep olnly columns with genes, start positions and chromosomes 
causal.pos.B.2 <- unique(causal.pos.B[,.(geneA, geneB, start.A, chr.A, chr.start, chr.id)])
colnames(causal.pos.B.2) <- c("geneA", "geneB", "start.A", "chr.A", "start.B", "chr.B")

## transform chromosome ids into numbers
# remove "chr" part of the chromosome name
causal.pos.B.2$chr.A <-gsub('chr', '', causal.pos.B.2$chr.A)
causal.pos.B.2$chr.B <-gsub('chr', '', causal.pos.B.2$chr.B)
colnames(causal.pos.B.2) <- c("geneA", "geneB", "start.A", "chr.A", "start.B", "chr.B")

# convert roman chromosome numbers to numbers
causal.pos.B.2$chr.A <- as.numeric(as.roman(causal.pos.B.2$chr.A))
causal.pos.B.2$chr.B <- as.numeric(as.roman(causal.pos.B.2$chr.B))

# order values
causal.pos.B.2.order <- causal.pos.B.2[order(chr.A, start.A, chr.B, start.B)]


# organize coordinates so that they are ordered by chromosome
# vector of chromosomes
vchr <- 1:16
# how much space will be separating chromosomes
separator <- 1e5

coordinates_plot <- sort_by_chr(vchr = vchr, causal.pos.B.2.order, separator = separator)
plot_sorted_coordinates(coordinates_plot, separator = separator)


#Create a function to generate a continuous color palette
rbPal <- colorRampPalette(c('blue','red'))


coordinates_plot_cor <- merge(coordinates_plot, find.effects[,.(geneA, geneB, cor)], by=c("geneA", "geneB"), all.x=T)
coordinates_plot_cor$cor <- abs(coordinates_plot_cor$cor)
coordinates_plot_cor <- coordinates_plot_cor[order(cor)]
coordinates_plot_cor$col <- rbPal(100)[as.numeric(cut(coordinates_plot_cor$cor,breaks = 100))]

plot_sorted_coordinates(coordinates_plot_cor, separator = separator, col = coordinates_plot_cor$col)
gradientLegend(format(round(coordinates_plot_cor$cor, 3), nsmall = 3), rbPal(100), inside=F, side=3)



coordinates_plot_cor.name <- unique(merge(coordinates_plot_cor, unique(genes_GO.table[,.(gene,symbol, gene.name, GO.term)]), 
      by.x="geneA", by.y="gene", all.x=T, allow.cartesian=TRUE))

genesA_names <- unique(coordinates_plot_cor.name[,.(geneA, start.A, chr.A, symbol, gene.name, GO.term)])
genesA_names[,.N, by=GO.term][order(-N)]

genesA_names2 <- unique(coordinates_plot_cor.name[,.(geneA, start.A, chr.A, symbol, gene.name, start.A)])
genesA_names2[,.N, by="gene.name"][order(-N)]
genesA_names2[order(-start.A)]
