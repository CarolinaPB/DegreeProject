library(data.table)
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))
genepos <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-20/gene_pos.gz")

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
# num chromosomes
nchr <- 16 

# how much space will be separating chromosomes
separator <- 62500
for (chr in 1:nchr){
  if (chr == 1){
    # coordinates for the causal genes
    final.coord.A <- data.table(chr, pos=causal.pos.B.2.order[chr.A == chr]$start.A)
    # coordinates for the affected genes
    final.coord.B <- data.table(chr, pos=causal.pos.B.2.order[chr.B == chr]$start.B)
    
  } else if (chr != 1) {
    previous <- chr-1
    # coordinates for the causal genes
    coord.A <- data.table(genechr, pos=max(final.coord.A[chr==previous]$pos) + causal.pos.B.2.order[chr.A == chr]$start.A + separator)
    # coordinates for the affected genes
    coord.B <- data.table(chr, pos=max(final.coord.B[chr==previous]$pos) + causal.pos.B.2.order[chr.B == chr]$start.B + separator)
    
    # add to the rest of the coordinates
    final.coord.A <- rbind(final.coord.A, coord.A)
    final.coord.B <- rbind(final.coord.B, coord.B)
  }
}

test <- data.table(causal.pos.B.2.order)

for (chr in 1:nchr){
  if (chr != 1) {
    previous <- chr-1
    test[chr.A == chr]$start.A <- max(test[chr.A==previous]$start.A) + test[chr.A == chr]$start.A + separator
    test[chr.B == chr]$start.B <- max(test[chr.B==previous]$start.B) + test[chr.B == chr]$start.B + separator
  }
}






# pdf("results/results_figures/affectedgenes_vs_causalgenes_position.pdf")
plot(test$start.A, test$start.B, pch=".", axes=F, xlab = "Causal gene position (chr)", ylab = "Affected gene position (chr)")

# add chromosome separators
for (ch in 1:nchr){
  abline(v= max(test[chr.A==ch]$start.A)+separator/2, col="lightblue", lty=2) 
  abline(h= max(test[chr.B==ch]$start.B)+separator/2, col="lightblue", lty=2) 
}

# add x and y axis chromosomes
axis(1, at=sapply(1:16, function(i){min(test[chr.A==i]$start.A) + (max(test[chr.A==i]$start.A) - min(test[chr.A==i]$start.A))/2}), labels=as.roman(1:16), tick=FALSE)
axis(2, at=sapply(1:16, function(i){min(test[chr.B==i]$start.B) + (max(test[chr.B==i]$start.B) - min(test[chr.B==i]$start.B))/2}), labels=as.roman(1:16), tick=FALSE)
# dev.off()
