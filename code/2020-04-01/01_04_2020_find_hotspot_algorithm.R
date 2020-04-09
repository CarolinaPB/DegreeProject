library(data.table)
load("results/coordinates_plot_cor.Rdata")
load("results/numlinks_from_gene.Rdata")
genepos <- fread("results/gene_pos.gz")

source("myfunctions.R")

coordinates_plot_links <- merge(coordinates_plot_cor, unique(numlinks_from_gene[,.(geneA, count.A)])[order(-count.A)], by="geneA")

causalgenes.pos.count1 <- unique(coordinates_plot_links[order(chr.A,start.A)][,.(geneA, start.A, end.A, chr.A, count.A)])

allgenes.location <- genepos

## transform chromosome ids into numbers
# remove "chr" part of the chromosome name
allgenes.location$chr.id <- gsub('chr', '', allgenes.location$chr.id)
# convert roman chromosome numbers to numbers
allgenes.location$chr.id <- as.numeric(as.roman(allgenes.location$chr.id))

# geneA with counts, sorted position and original position
causalgenes.pos.count <- merge(causalgenes.pos.count1, allgenes.location, by.x=c("geneA", "chr.A"), by.y=c("gene", "chr.id"))
causalgenes.pos.count <- causalgenes.pos.count[order(chr.A, start.A, end.A)]

# For each chromosome get the sequence of values that follow a certain condition
# condition == count.A>10
# True - causal genes affect 10 or more genes
# False - causal genes affect less than 10 genes

# create list of tables - one for each chr
# each table has the runs for the condition count.A>10

rle.res.list <- rle_causalgenes(causalgenes.pos.count, lim = 10)


# find hotspots and plot them
# pdf(file = "results/figures/my_causal_hotspots.pdf", onefile = T)
plot_sorted_coordinates(coordinates_plot_cor, separator = separator, col = coordinates_plot_cor[chr.A==chr]$col)

hot <- find_hotspots(rle.list = rle.res.list, coordinates_plot = coordinates_plot_cor, causalgenes = causalgenes.pos.count, lim = 10, plt = T)
# dev.off()

# save the hotspot intervals in a grange object (also available in a table)
granges_myhotspots <- GRanges(seqnames = hot$chr,
                            ranges=IRanges(start=hot$start,
                                           end=hot$end))
