library(data.table)
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))
genepos <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-20/gene_pos.gz")

source(file = "code/myfunctions.R")
# Add chomosome and start position to each gene
causal.pos.A <- merge(find.effects_TF, genepos, by.x="geneA", by.y="gene", all.x=T)
colnames(causal.pos.A) <- c("geneA", "geneB", "eqtl.A", "eqtl.B", "A->B", "B->A", "strand.A", "start.A","end.A", "chr.A")
causal.pos.B <- merge(causal.pos.A, genepos, by.x="geneB", by.y="gene", all.x=T)

# Keep olnly columns with genes, start positions and chromosomes 
causal.pos.B.2 <- unique(causal.pos.B[,.(geneA, geneB, start.A, chr.A, start, chr)])

## transform chromosome ids into numbers
# remove "chr" part of the chromosome name
causal.pos.B.2$chr.A <-gsub('chr', '', causal.pos.B.2$chr.A)
causal.pos.B.2$chr <-gsub('chr', '', causal.pos.B.2$chr)
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

