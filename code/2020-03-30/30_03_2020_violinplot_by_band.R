plot_sorted_coordinates(coordinates_plot_cor, separator = separator, col = coordinates_plot_cor[chr.A==chr]$col)


# plot chr 7
chr <- 7
## plot
plot_sorted_coordinates(coordinates_plot_cor[chr.A==chr], separator = separator, col = coordinates_plot_cor[chr.A==chr]$col)

## define the bands
chr7.causal.pos.counts <- unique(coordinates_plot_links[chr.A==chr][order(start.A)][,.(geneA, start.A, end.A, chr.A, count.A)])
### first band
chr7.band1.min <- 4152338 # start pos of the first gene in the band
chr7.band1.max <- 4297719 # end pos of the last gene in the band
abline(v=c(chr7.band1.min, chr7.band1.max), col="lightblue")

### get genomic ranges for chr 7
granges.chr7 <- GRanges(seqnames = chr7.causal.pos.counts$chr.A,
                     ranges=IRanges(start=chr7.causal.pos.counts$start.A,
                                    end = chr7.causal.pos.counts$end.A))
names(granges.chr7) <- chr7.causal.pos.counts$geneA
granges.chr7 <- unique(granges.chr7)

### get genomic ranges limits for the band
gr.chr7.band1.lims <- GRanges(seqnames = chr,
                        ranges=IRanges(start=chr7.band1.min,
                                       end = chr7.band1.max))

### keep only the genes inside the band
chr7.genes.band1 <- subsetByOverlaps(granges.chr7, gr.chr7.band1.lims)
chr7.genes.band1.names <- names(chr7.genes.band1)

## Violin plot reads genes in first band
genes_to_plot.chr7.band1 <- data.table(genes = chr7.genes.band1.names, chr)

### get read counts for those genes
chr7.band1.reads <- get_readcounts_toplot(genes_to_plot = genes_to_plot.chr7.band1, chr=7, causal = "causal")

### plot
p.sub <- ggplot(chr7.band1.reads, aes(x=gene, y=counts, fill=as.factor(chr))) + 
  geom_violin() + 
  labs(title=paste0("Read counts chr",chr))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

## second band
chr7.band2.min <- 4455928
chr7.band2.max <- 4608279
abline(v=c(chr7.band2.min, chr7.band2.max), col="lightblue")

gr.chr7.band2.lims <- GRanges(seqnames = chr,
                              ranges=IRanges(start=chr7.band2.min,
                                             end = chr7.band2.max))

### keep only the genes inside the band
chr7.genes.band2 <- subsetByOverlaps(granges.chr7, gr.chr7.band2.lims)
chr7.genes.band2.names <- names(chr7.genes.band2)

## Violin plot reads genes in second band
genes_to_plot.chr7.band2 <- data.table(genes = chr7.genes.band2.names, chr)

### get read counts for those genes
chr7.band2.reads <- get_readcounts_toplot(genes_to_plot = genes_to_plot.chr7.band2, chr=7, causal = "causal")

### plot
p.sub <- ggplot(chr7.band1.reads, aes(x=gene, y=counts, fill=as.factor(chr))) + 
  geom_violin() + 
  labs(title=paste0("Read counts chr",chr))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


## Plot reads for genes in both bands of chr7
genes_to_plot.chr7.all_bands <- rbind(data.table(genes_to_plot.chr7.band1, band="1"), data.table(genes_to_plot.chr7.band2,  band="2"))
chr7.allbands.reads <- get_readcounts_toplot(genes_to_plot = genes_to_plot.chr7.all_bands, chr=7, causal = "causal")

p.sub <- ggplot(chr7.allbands.reads, aes(x=gene, y=counts, fill=band)) + 
  geom_violin() + 
  labs(title=paste0("Read counts chr",chr))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 



######
allgenes.location <- genepos

## transform chromosome ids into numbers
# remove "chr" part of the chromosome name
allgenes.location$chr.id <- gsub('chr', '', allgenes.location$chr.id)
# convert roman chromosome numbers to numbers
allgenes.location$chr.id <- as.numeric(as.roman(allgenes.location$chr.id))

# genes in chr7 that have no causal effect
chr7_notcausal <- allgenes.location[chr.id==7 & !gene %in% genes_to_plot.chr7.all_bands$genes][order(chr.start)]

# see if they are inside bands
## define bands
coordinates_plot_oldpos <- merge(coordinates_plot_links, allgenes.location[,.(gene, chr.start, chr.end)], by.y="gene", by.x="geneA")

chr7.causal.pos.counts <- unique(coordinates_plot_oldpos[chr.A==chr][order(start.A)][,.(geneA, start.A, end.A, chr.A, count.A,chr.start, chr.end)])

# gene YGL236C
chr7.band1.min <- 53787 # start pos of the first gene in the band
# gene YGL166W
chr7.band1.max <- 191806 # end pos of the last gene in the band

### get genomic ranges for chr 7 (original positions, not the sorted positions)
granges.chr7 <- GRanges(seqnames = chr7.causal.pos.counts$chr.A,
                        ranges=IRanges(start=chr7.causal.pos.counts$chr.start,
                                       end = chr7.causal.pos.counts$chr.end))
names(granges.chr7) <- chr7.causal.pos.counts$geneA
granges.chr7 <- unique(granges.chr7)

### get genomic ranges limits for the band
gr.chr7.band1.lims <- GRanges(seqnames = chr,
                              ranges=IRanges(start=chr7.band1.min,
                                             end = chr7.band1.max))

### keep only the genes inside the band
chr7.genes.band1 <- subsetByOverlaps(granges.chr7, gr.chr7.band1.lims)
chr7.genes.band1.names <- names(chr7.genes.band1)

## Violin plot reads genes in first band
genes_to_plot.chr7.band1 <- data.table(genes = chr7.genes.band1.names, chr, band=1)

### get read counts for those genes
chr7.band1.reads <- get_readcounts_toplot(genes_to_plot = genes_to_plot.chr7.band1, chr=7, causal = "causal")

### plot
p.sub <- ggplot(chr7.band1.reads, aes(x=gene, y=counts, fill=as.factor(chr))) + 
  geom_violin() + 
  labs(title=paste0("Read counts chr",chr))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

## second band
# gene YGL081W
chr7.band2.min <- 357377
#YGR003W
chr7.band2.max <- 502366

gr.chr7.band2.lims <- GRanges(seqnames = chr,
                              ranges=IRanges(start=chr7.band2.min,
                                             end = chr7.band2.max))

### keep only the genes inside the band
chr7.genes.band2 <- subsetByOverlaps(granges.chr7, gr.chr7.band2.lims)
chr7.genes.band2.names <- names(chr7.genes.band2)

## Violin plot reads genes in second band
genes_to_plot.chr7.band2 <- data.table(genes = chr7.genes.band2.names, chr, band=2)

### get read counts for those genes
chr7.band2.reads <- get_readcounts_toplot(genes_to_plot = genes_to_plot.chr7.band2, chr=7, causal = "causal")

### plot
p.sub <- ggplot(chr7.band2.reads, aes(x=gene, y=counts, fill=as.factor(chr))) + 
  geom_violin() + 
  labs(title=paste0("Read counts chr",chr))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


## Plot reads for genes in both bands of chr7
genes_to_plot.chr7.all_bands <- rbind(data.table(genes_to_plot.chr7.band1), data.table(genes_to_plot.chr7.band2))
chr7.allbands.reads <- get_readcounts_toplot(genes_to_plot = genes_to_plot.chr7.all_bands, chr=7, causal = "causal")

p.sub <- ggplot(chr7.allbands.reads, aes(x=gene, y=counts, fill=as.factor(band))) + 
  geom_violin() + 
  labs(title=paste0("Read counts chr",chr))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

### Genomic ranges for genes with no causal effects on chr 7
granges.chr7.notcausal <- GRanges(seqnames = chr7_notcausal$chr.id,
                        ranges=IRanges(start=chr7_notcausal$chr.start,
                                       end = chr7_notcausal$chr.end))
names(granges.chr7.notcausal) <- chr7_notcausal$gene
granges.chr7.notcausal <- unique(granges.chr7.notcausal)

# genes with no causal effect on band 1 of chr7
chr7.notcausal.band1 <- subsetByOverlaps(granges.chr7.notcausal, gr.chr7.band1.lims)
chr7.notcausal.band1.names <- names(chr7.notcausal.band1)
# genes with no causal effect on band 2 of chr7
chr7.notcausal.band2 <- subsetByOverlaps(granges.chr7.notcausal, gr.chr7.band2.lims)
chr7.notcausal.band2.names <- names(chr7.notcausal.band2)

chr7.notcausal.band1 <- data.table(genes = chr7.notcausal.band1.names, chr, band=1)
chr7.notcausal.band2 <- data.table(genes = chr7.notcausal.band2.names, chr, band=2)

notcausal.chr7.all_bands <- rbind(chr7.notcausal.band1, genes_to_plot.chr7.band2)
notcausal.chr7.allbands.reads <- get_readcounts_toplot(genes_to_plot = notcausal.chr7.all_bands, chr=7, causal = "no")

p.sub <- ggplot(notcausal.chr7.allbands.reads, aes(x=gene, y=counts, fill=as.factor(band))) + 
  geom_violin() + 
  labs(title=paste0("Read counts chr",chr, " - not causal genes"))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

chr7.violin <- rbind(chr7.allbands.reads, notcausal.chr7.allbands.reads)

p.sub <- ggplot(chr7.violin[gene %in% unique(chr7.violin$gene)[1:10]], aes(x=gene, y=counts, fill=as.factor(causal))) + 
  geom_violin() + 
  ylim(0,2000)+
  labs(title=paste0("Read counts chr",chr))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

p.sub <- ggplot(chr7.violin, aes(x=gene, y=counts, fill=as.factor(causal))) + 
  geom_violin() + 
  labs(title=paste0("Read counts chr",chr))
p.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1)) 


####
# pick random genes that are not in any hotspot
## define all hotspots bands
causal.hotspots <- list()
causal.hotspots <- gr.chr7.band1.lims
causal.hotspots <- c(causal.hotspots, gr.chr7.band2.lims)
names(causal.hotspots)[1] <- "chr7.band1"
names(causal.hotspots)[2] <- "chr7.band2"


# chr8
chr <- 8

chr8.causal.pos.counts <- unique(coordinates_plot_oldpos[chr.A==chr][order(start.A)][,.(geneA, start.A, end.A, chr.A, count.A, chr.start, chr.end)])
plot_sorted_coordinates(coordinates_plot_cor[chr.A==chr], 
                        separator = separator, 
                        col = coordinates_plot_cor[chr.A==chr]$col)

### first band
#YHL027W
chr8.band1.min <- 51111 # start pos of the first gene in the band
# YHR064C
chr8.band1.max <- 227141 # end pos of the last gene in the band

### get genomic ranges for chr 8
granges.chr8 <- GRanges(seqnames = chr8.causal.pos.counts$chr.A,
                        ranges=IRanges(start=chr8.causal.pos.counts$start.A,
                                       end = chr8.causal.pos.counts$end.A))
names(granges.chr8) <- chr8.causal.pos.counts$geneA
granges.chr8 <- unique(granges.chr8)

### get genomic ranges limits for the band
gr.chr8.band1.lims <- GRanges(seqnames = chr,
                              ranges=IRanges(start=chr8.band1.min,
                                             end = chr8.band1.max))

causal.hotspots <- c(causal.hotspots, gr.chr8.band1.lims)
names(causal.hotspots)[3] <- "chr8.band1"

# chr12
chr <- 12

chr12.causal.pos.counts <- unique(coordinates_plot_oldpos[chr.A==chr][order(start.A)][,.(geneA, start.A, end.A, chr.A, count.A, chr.start, chr.end)])
plot_sorted_coordinates(coordinates_plot_cor[chr.A==chr], 
                        separator = separator, 
                        col = coordinates_plot_cor[chr.A==chr]$col)

### get genomic ranges for chr 12
granges.chr12 <- GRanges(seqnames = chr12.causal.pos.counts$chr.A,
                         ranges=IRanges(start=chr12.causal.pos.counts$start.A,
                                        end = chr12.causal.pos.counts$end.A))
names(granges.chr12) <- chr12.causal.pos.counts$geneA
granges.chr12 <- unique(granges.chr12)

### first band
# YLL023C
chr12.band1.min <- 97997 # start pos of the first gene in the band
# YLR001C
chr12.band1.max <- 153977 # end pos of the last gene in the band

### get genomic ranges limits for the band
gr.chr12.band1.lims <- GRanges(seqnames = chr,
                              ranges=IRanges(start=chr12.band1.min,
                                             end = chr12.band1.max))

### Second band
# YLR215C
chr12.band2.min <- 570776 # start pos of the first gene in the band
# YLR351C
chr12.band2.max <- 830364 # end pos of the last gene in the band

### get genomic ranges limits for the band
gr.chr12.band2.lims <- GRanges(seqnames = chr,
                               ranges=IRanges(start=chr12.band2.min,
                                              end = chr12.band2.max))

### Third band
# YLR371W
chr12.band3.min <- 862714 # start pos of the first gene in the band
# YLR438W
chr12.band3.max <- 1013775 # end pos of the last gene in the band

### get genomic ranges limits for the band
gr.chr12.band3.lims <- GRanges(seqnames = chr,
                               ranges=IRanges(start=chr12.band3.min,
                                              end = chr12.band3.max))


causal.hotspots <- c(causal.hotspots, gr.chr12.band1.lims)
names(causal.hotspots)[4] <- "chr12.band1"
causal.hotspots <- c(causal.hotspots, gr.chr12.band2.lims)
names(causal.hotspots)[5] <- "chr12.band2"
causal.hotspots <- c(causal.hotspots, gr.chr12.band3.lims)
names(causal.hotspots)[6] <- "chr12.band3"

# chr13
chr <- 13

chr13.causal.pos.counts <- unique(coordinates_plot_oldpos[chr.A==chr][order(start.A)][,.(geneA, start.A, end.A, chr.A, count.A, chr.start, chr.end)])
plot_sorted_coordinates(coordinates_plot_cor[chr.A==chr], 
                        separator = separator, 
                        col = coordinates_plot_cor[chr.A==chr]$col)

### get genomic ranges for chr 13
granges.chr13 <- GRanges(seqnames = chr13.causal.pos.counts$chr.A,
                         ranges=IRanges(start=chr13.causal.pos.counts$start.A,
                                        end = chr13.causal.pos.counts$end.A))
names(granges.chr13) <- chr13.causal.pos.counts$geneA
granges.chr13 <- unique(granges.chr13)

### first band
# YML020W
chr13.band1.min <- 97997 # start pos of the first gene in the band
# YMR068W
chr13.band1.max <- 407584 # end pos of the last gene in the band

### get genomic ranges limits for the band
gr.chr13.band1.lims <- GRanges(seqnames = chr,
                               ranges=IRanges(start=chr13.band1.min,
                                              end = chr13.band1.max))

### Second band
# YLR215C
chr13.band2.min <- 570776 # start pos of the first gene in the band
# YLR351C
chr13.band2.max <- 830364 # end pos of the last gene in the band

### get genomic ranges limits for the band
gr.chr13.band2.lims <- GRanges(seqnames = chr,
                               ranges=IRanges(start=chr13.band2.min,
                                              end = chr13.band2.max))

### Third band
# YLR371W
chr13.band3.min <- 862714 # start pos of the first gene in the band
# YLR438W
chr13.band3.max <- 1013775 # end pos of the last gene in the band

### get genomic ranges limits for the band
gr.chr13.band3.lims <- GRanges(seqnames = chr,
                               ranges=IRanges(start=chr13.band3.min,
                                              end = chr13.band3.max))


causal.hotspots <- c(causal.hotspots, gr.chr13.band1.lims)
names(causal.hotspots)[7] <- "chr13.band1"
causal.hotspots <- c(causal.hotspots, gr.chr13.band2.lims)
names(causal.hotspots)[8] <- "chr13.band2"
causal.hotspots <- c(causal.hotspots, gr.chr13.band3.lims)
names(causal.hotspots)[9] <- "chr13.band3"

# chr14
chr <- 14

chr14.causal.pos.counts <- unique(coordinates_plot_oldpos[chr.A==chr][order(start.A)][,.(geneA, start.A, end.A, chr.A, count.A, chr.start, chr.end)])
plot_sorted_coordinates(coordinates_plot_cor[chr.A==chr], 
                        separator = separator, 
                        col = coordinates_plot_cor[chr.A==chr]$col)

### get genomic ranges for chr 14
granges.chr14 <- GRanges(seqnames = chr14.causal.pos.counts$chr.A,
                         ranges=IRanges(start=chr14.causal.pos.counts$start.A,
                                        end = chr14.causal.pos.counts$end.A))
names(granges.chr14) <- chr14.causal.pos.counts$geneA
granges.chr14 <- unique(granges.chr14)

### first band
# YML020W
chr14.band1.min <- 97997 # start pos of the first gene in the band
# YNL010W
chr14.band1.max <- 614360 # end pos of the last gene in the band

### get genomic ranges limits for the band
gr.chr14.band1.lims <- GRanges(seqnames = chr,
                               ranges=IRanges(start=chr14.band1.min,
                                              end = chr14.band1.max))


causal.hotspots <- c(causal.hotspots, gr.chr14.band1.lims)
names(causal.hotspots)[10] <- "chr14.band1"

# chr15
chr <- 15

chr15.causal.pos.counts <- unique(coordinates_plot_oldpos[chr.A==chr][order(start.A)][,.(geneA, start.A, end.A, chr.A, count.A, chr.start, chr.end)])
plot_sorted_coordinates(coordinates_plot_cor[chr.A==chr], 
                        separator = separator, 
                        col = coordinates_plot_cor[chr.A==chr]$col)

### get genomic ranges for chr 15
granges.chr15 <- GRanges(seqnames = chr15.causal.pos.counts$chr.A,
                         ranges=IRanges(start=chr15.causal.pos.counts$start.A,
                                        end = chr15.causal.pos.counts$end.A))
names(granges.chr15) <- chr15.causal.pos.counts$geneA
granges.chr15 <- unique(granges.chr15)

### first band
# YOL154W
chr15.band1.min <- 34658 # start pos of the first gene in the band
# YOL049W
chr15.band1.max <- 249534 # end pos of the last gene in the band

### get genomic ranges limits for the band
gr.chr15.band1.lims <- GRanges(seqnames = chr,
                               ranges=IRanges(start=chr15.band1.min,
                                              end = chr15.band1.max))


causal.hotspots <- c(causal.hotspots, gr.chr15.band1.lims)
names(causal.hotspots)[11] <- "chr15.band1"


# Get Granges for all genes
allgenes.location[is.na(chr.id)]$chr.id <- 20
granges_allgenes <- GRanges(seqnames = allgenes.location$chr.id,
                            ranges=IRanges(start=allgenes.location$chr.start,
                                           end = allgenes.location$chr.end))

# get genes that are not in hotspots
genes_in_hotspots <- subsetByOverlaps(granges_allgenes, causal.hotspots)
genes_notin_hotspots <- allgenes.location[!gene %in%  names(genes_in_hotspots)]$gene

v <- chr12.causal.pos.counts$count.A > 10
r <- rle(v)
