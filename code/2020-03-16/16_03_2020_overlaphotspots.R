library(GenomicRanges)

# Add chomosome and start position to each gene
causal.pos.A <- merge(find.effects_TF, genepos, by.x="geneA", by.y="gene", all.x=T)
colnames(causal.pos.A) <- c("geneA", "geneB", "eqtl.A", "eqtl.B", "A->B", "B->A", "strand.A", "start.A","end.A", "chr.A")
causal.pos.B <- merge(causal.pos.A, genepos, by.x="geneB", by.y="gene", all.x=T)


# Keep olnly columns with genes, start positions and chromosomes
gene.location <- unique(causal.pos.B[,.(geneA, geneB, start.A,end.A, chr.A, chr.start,chr.end, chr.id, strand.A, chr.strand)])
colnames(gene.location) <- c("geneA", "geneB", "start.A","end.A", "chr.A", "start.B","end.B", "chr.B", "strand.A", "strand.B")

## transform chromosome ids into numbers
# remove "chr" part of the chromosome name
gene.location$chr.A <-gsub('chr', '', gene.location$chr.A)
gene.location$chr.B <-gsub('chr', '', gene.location$chr.B)

# convert roman chromosome numbers to numbers
gene.location$chr.A <- as.numeric(as.roman(gene.location$chr.A))
gene.location$chr.B <- as.numeric(as.roman(gene.location$chr.B))

# order values
gene.location.order <- gene.location[order(chr.A, start.A, chr.B, start.B)]


granges.A <- GRanges(seqnames = gene.location.order$chr.A,
        ranges=IRanges(start=gene.location.order$start.A,
                       end = gene.location.order$end.A),
        strand= gene.location.order$strand.A)
names(granges.A) <- gene.location.order$geneA
granges.A <- unique(granges.A)

granges.B <- GRanges(seqnames = gene.location.order$chr.B,
                     ranges=IRanges(start=gene.location.order$start.B,
                                    end = gene.location.order$end.B),
                     strand= gene.location.order$strand.B)
names(granges.B) <- gene.location.order$geneB
granges.B <- unique(granges.B)


granges.hotspots <- GRanges(seqnames = hotspot_data.interval$chromosome,
                            ranges=IRanges(start=hotspot_data.interval$bootstrapIntervalLeft,
                                           end=hotspot_data.interval$bootstrapIntervalRight))
names(granges.hotspots) <- hotspot_data.interval$hotspotMarker


# genes that are in hotspots
overlapping_genesA <- subsetByOverlaps(granges.A, granges.hotspots)
genes_in_hotspot <- unique(gene.location.order[geneA %in% names(overlapping_genesA)][,.(geneA, start.A, end.A, chr.A, strand.A)])


# split the Grange object by chromosome
splitchr <- split(granges.A, seqnames(granges.A))

findOverlaps(granges.A, granges.hotspots, ignore.strand=T)
# data.table(countOverlaps(granges.A, granges.hotspots))
subsetByOverlaps(granges.A, granges.hotspots)
countOverlaps(granges.A, granges.hotspots)
# how many genes are in each hotspot?
coverage(granges.A)
max(coverage(granges.A))
mean(coverage(granges.A))

ranges(granges.A)
seqnames(granges.A)
start(granges.A)
end(granges.A)

over <- findOverlaps(granges.A, granges.hotspots, ignore.strand=T)
queryHits(over)
subjectHits(over)

genes_in_hotspot$hotspot <- names(granges.hotspots[subjectHits(over)])



