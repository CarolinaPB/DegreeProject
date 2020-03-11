sortmap <- function (chrom, map, delta = 1) #Stolen from GenABEL
{
  chnum <- as.numeric(as.factor(chrom))
  ix <- order(chnum, map)
  map <- map[ix]
  off <- c(0, map[1:(length(map) - 1)])
  off <- map - off
  off[which(off <= 0)] <- delta
  cummap <- cumsum(off)
  out <- list()
  out$ix <- ix
  out$cummap <- cummap
  out$chnum <- chnum
  out
}

genepos <- fread("results/gene_pos.gz")
if (!exists("find.effects_TF")){
  find.effects_TF <- fread("results/findeffects_TF_newparams.gz")
}

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

datA <- data.frame(unique(causal.pos.B.2[,.(geneA, chr.A, start.A)]))
datB <- data.frame(unique(causal.pos.B.2[,.(geneB, chr.B, start.B)]))
# dat <- data.frame(chr = sample(c('chr1', 'chr2', 'chr3'), size = 100, replace = T), 
#                   pos = sample(1:1e6, size = 100))
dat.mapA <- sortmap(chrom = datA$chr.A, map = datA$start.A)
dat.mapB <- sortmap(chrom = datB$chr.B, map = datB$start.B)
# dat.map$ix     #This are the indices that sorts dat by chr and position
# dat.map$cummap #The cummulative positions. This refers to the positions re-ordered by dat.map$ix
#Example
AA <- cbind(datA[dat.mapA$ix, ], dat.mapA$cummap)
BB <- cbind(datB[dat.mapB$ix, ], dat.mapB$cummap)

withcumsum1 <- merge(causal.pos.B.2, AA, by=c("geneA", "start.A", "chr.A"))
withcumsum <- merge(withcumsum1, BB, by=c("geneB", "start.B", "chr.B"))

plot(withcumsum[chr.A == 12 &chr.B == 3]$`dat.mapA$cummap`, withcumsum[chr.A == 12 &chr.B == 3]$`dat.mapB$cummap`, pch=".")

# get the same thing with sortmap and with my sorting function 
# genesA_start_order <- unique(coordinates_plot_cor[chr.A == 12 &chr.B == 3][order(start.A)][, .(geneA, start.A)])
# 
# plot_sorted_coordinates(
#   coordinates_plot_cor[chr.A == 12&chr.B == 3],
#   separator = separator,
#   col = coordinates_plot_cor$col,
# )
# abline(v=c(7463067-1500, 7526178+1500), col="pink")