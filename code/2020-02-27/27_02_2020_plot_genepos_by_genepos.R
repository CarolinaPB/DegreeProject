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

respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))
genepos <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-20/gene_pos.gz")

causal.pos.A <- merge(find.effects_TF, genepos, by.x="geneA", by.y="gene", all.x=T)
colnames(causal.pos.A) <- c("geneA", "geneB", "eqtl.A", "eqtl.B", "A->B", "B->A", "strand.A", "start.A","end.A", "chr.A")
causal.pos.B <- merge(causal.pos.A, genepos, by.x="geneB", by.y="gene", all.x=T)
plot(causal.pos.B$start.A, causal.pos.B$start, pch=".")

causal.pos.A.2 <- unique(causal.pos.A[,.(geneA, start.A, chr.A)])
causal.pos.B.2 <- unique(causal.pos.B[,.(geneA, geneB, start.A, chr.A, start, chr)])


causal.pos.B.2$chr.A <-gsub('chr', '', causal.pos.B.2$chr.A)
causal.pos.B.2$chr <-gsub('chr', '', causal.pos.B.2$chr)
colnames(causal.pos.B.2) <- c("geneA", "geneB", "start.A", "chr.A", "start.B", "chr.B")

# convert roman chromosome numbers to numbers
causal.pos.B.2$chr.A <- as.numeric(as.roman(causal.pos.B.2$chr.A))
causal.pos.B.2$chr.B <- as.numeric(as.roman(causal.pos.B.2$chr.B))

causal.pos.B.2.order <- causal.pos.B.2[order(chr.A, start.A, chr.B, start.B)]

max.startchrA <- causal.pos.B.2.order[, .(max=max(start.A)), by="chr.A"]
max.startchrB <- causal.pos.B.2.order[, .(max=max(start.B)), by="chr.B"]
maxpos_chr <- merge(max.startchrA, max.startchrB, by.x="chr.A", by.y = "chr.B")
colnames(maxpos_chr) <- c("chr", "causal", "affected")

nchr <- 16
final.coord.A <- data.table(chr=numeric(), pos=numeric())
final.coord.B <- data.table(chr=numeric(), pos=numeric())
separator <- 500000
for (chr in 1:nchr){
  if (chr == 1){
    final.coord.A <- data.table(chr, pos=causal.pos.B.2.order[chr.A == chr]$start.A)
    final.coord.B <- data.table(chr, pos=causal.pos.B.2.order[chr.B == chr]$start.B)
    
  } else if (chr != 1) {
    previous <- chr-1
    coord.A <- data.table(chr, pos=max(final.coord.A[chr==previous]$pos) + causal.pos.B.2.order[chr.A == chr]$start.A + separator)
    coord.B <- data.table(chr, pos=max(final.coord.B[chr==previous]$pos) + causal.pos.B.2.order[chr.B == chr]$start.B + separator)
    
    
    final.coord.A <- rbind(final.coord.A, coord.A)
    final.coord.B <- rbind(final.coord.B, coord.B)
  }
}

plot(final.coord.A$pos, final.coord.B$pos, pch=".", axes=F, xlab = "Causal gene position (chr)", ylab = "Affected gene position (chr)")

for (ch in 1:nchr){
  abline(v= max(final.coord.A[chr==ch]$pos)+separator/2, col="lightblue", lty=2) 
  abline(h= max(final.coord.B[chr==ch]$pos)+separator/2, col="lightblue", lty=2) 
}



axis(1, at=sapply(1:16, function(i){min(final.coord.A[chr==i]$pos) + (max(final.coord.A[chr==i]$pos) - min(final.coord.A[chr==i]$pos))/2}), labels=as.roman(1:16), tick=FALSE)
axis(2, at=sapply(1:16, function(i){min(final.coord.B[chr==i]$pos) + (max(final.coord.B[chr==i]$pos) - min(final.coord.B[chr==i]$pos))/2}), labels=as.roman(1:16), tick=FALSE)


# plot(causal.pos.B.2.order$start.A, causal.pos.B.2.order$start.B, pch = ".")
# abline(h = max(causal.pos.B.2.order[chr.A ==1]$start.A))
# abline(v = max(causal.pos.B.2.order[chr.A ==1]$start.A))
# abline(h = max(causal.pos.B.2.order[chr.A ==1]$start.B))
# abline(v = max(causal.pos.B.2.order[chr.A ==1]$start.B))
# 
# plot(causal.pos.B.2.order[chr.A ==1]$start.A, pch = ".")




testA <- sortmap(as.numeric(as.roman(causal.pos.B.2$chr.A)), map = causal.pos.B.2$start.A)
testB <- sortmap(as.numeric(as.roman(causal.pos.B.2$chr.B)), map = causal.pos.B.2$start.B)
plot(testA$ix, testB$ix, pch=".", xlab="causal", ylab="affected", xaxt="n", yaxt="n", xaxs="i", yaxs="i")

chromosomeDividers <- c(0, 10000, 209653, 1360022, 2891955, 3468829, 3738990, 4829930, 5392573, 5832461, 6578212, 7245028, 8323205, 9247636, 10031969, 11123260, 12071326)


for(i in 1:length(chromosomeDividers)){
  abline(v=chromosomeDividers[i], col="red", lty=1)
  abline(h=chromosomeDividers[i], col="red", lty=1)
}
axis(1, at=chromosomeDividers, tick=FALSE)
axis(2, at=sapply(1:16, function(i){chromosomeDividers[i] + ((chromosomeDividers[i+1] - chromosomeDividers[i])/2)}), labels=as.roman(1:16), tick=FALSE)


chrom <- causal.pos.B$chr.A
map <- causal.pos.B$start.A

plot(causal.pos.B.2[as.numeric(as.roman(chr.A))==1 & as.numeric(as.roman(chr.B))==1][
  order(as.numeric(as.roman(chr.A)))]$start.A, causal.pos.B.2[as.numeric(as.roman(chr.A))==1 & 
                                                              as.numeric(as.roman(chr.B))==1][order(as.numeric(as.roman(chr.A)))]$start.B, pch=".")


