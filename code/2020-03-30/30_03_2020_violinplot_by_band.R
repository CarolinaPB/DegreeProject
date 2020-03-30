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

## Violin plot reads genes in first band
genes_to_plot.chr7.band2 <- data.table(genes = chr7.genes.band2.names, chr)

### get read counts for those genes
chr7.band1.reads <- get_readcounts_toplot(genes_to_plot = genes_to_plot.chr7.band2, chr=7, causal = "causal")

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



