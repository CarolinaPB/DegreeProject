library(data.table)
library("tidyr")

genepos <- fread("results/gene_pos.gz")
find.effects_TF <- fread("results/findeffects_TF_newparams.gz")

# Add chomosome and start position to each gene
causal.pos.A <- merge(find.effects_TF, genepos, by.x="geneA", by.y="gene", all.x=T)
colnames(causal.pos.A) <- c("geneA", "geneB", "eqtl.A", "eqtl.B", "A->B", "B->A", "strand.A", "start.A","end.A", "chr.A")
causal.pos.B <- merge(causal.pos.A, genepos, by.x="geneB", by.y="gene", all.x=T)

# separate eqtl name into chromosome number and position of the eqtl
causal.pos.eqtlA <- separate(causal.pos.B, eqtl.A, c("eqtlA.chr", "eqtlA.pos", "eqtlA.other"), sep = "[:_]+", remove = F,
                             convert = T, extra = "warn", fill = "warn")
causal.pos.eqtlB <- separate(causal.pos.eqtlA, eqtl.B, c("eqtlB.chr", "eqtlB.pos", "eqtlB.other"), sep = "[:_]+", remove = F,
                             convert = T, extra = "warn", fill = "warn")

# tidy up
colstoremove <-  c("eqtlA.other", "eqtlB.other", "strand.A", "chr.strand")
causal.pos.eqtlB[,paste0(colstoremove):=NULL]
setnames(causal.pos.eqtlB, c("geneB", "geneA","eqtl.A","eqtlA.chr", "eqtlA.pos", 
                             "eqtl.B","eqtlB.chr", "eqtlB.pos", "A->B", "B->A", 
                             "start.A","end.A", "chr.A", "start.B","end.B","chr.B"))
setcolorder(causal.pos.eqtlB, c("geneA", "geneB", "eqtl.A", "eqtl.B","chr.A", 
                                "start.A","end.A","chr.B", "start.B","end.B","eqtlA.chr", 
                                "eqtlA.pos", "eqtlB.chr", "eqtlB.pos"))

# remove "chr" from chromosome number
cols = c("chr.A","chr.B", "eqtlA.chr", "eqtlB.chr")   # define which columns to work with
causal.pos.eqtlB[ , (cols) := lapply(.SD, function(x) {gsub("chr", "", x)}), .SDcols = cols] # replace
# "translate" into roman numerals into numeric
causal.pos.eqtlB[ , (cols) := lapply(.SD, function(x) {as.numeric(as.roman(x))}), .SDcols = cols]

# calculate distance between gene and eqtl
causal.pos.eqtlB[,("dist.A"):=start.A-eqtlA.pos]
causal.pos.eqtlB[,("dist.B"):=start.B-eqtlB.pos]

# plot distance between gene and eqtl
plot(unique(causal.pos.eqtlB[,.(geneB, eqtl.B, dist.B)])[order(-dist.B)]$dist.B, pch=".", col="red",
     cex=2, ylab="distance", xlab="gene-eqtl pair", main="Distance between eqtl and gene")
points(unique(causal.pos.eqtlB[,.(geneA, eqtl.A, dist.A)])[order(-dist.A)]$dist.A, pch=".", col="blue", cex=2)


# Get genes that have eqtls "inside" the gene
unique(causal.pos.eqtlB[,.(geneA, eqtl.A, dist.A)])[order(-dist.A)][abs(dist.A) <= 100]

# negative distances mean that the eqtl is "inside" or after the gene
# positive numbers mean that the eqtl is before the gene 
# genesA that have an eqtl "inside"
distance.A <- unique(causal.pos.eqtlB[,.(geneA, eqtl.A, chr.A, start.A, end.A, eqtlA.chr,eqtlA.pos, dist.A)])
distance.B <- unique(causal.pos.eqtlB[,.(geneB, eqtl.B, chr.B, start.B, end.B, eqtlB.chr,eqtlB.pos, dist.B)])

# genes that have an eqtl "inside" the gene 
distA.inside <- distance.A[start.A < eqtlA.pos & eqtlA.pos < end.A][order(-dist.A)]
# genesA that have an eqtl before their start position
distA.before <- distance.A[start.A > eqtlA.pos][order(-dist.A)]
# genes that have an eqtl after their end position
distA.after <- distance.A[eqtlA.pos > end.A][order(-dist.A)]


# genes that have an eqtl "inside" the gene 
distB.inside <- distance.B[start.B < eqtlB.pos & eqtlB.pos < end.B][order(-dist.B)]
# genesA that have an eqtl before their start position
distB.before <- distance.B[start.B > eqtlB.pos][order(-dist.B)]
# genes that have an eqtl after their end position
distB.after <- distance.B[eqtlB.pos > end.B][order(-dist.B)]

# Plot number of genes that have eqtls before the gene start site, within the gene or after the gene's end site
toplotA <- c(before=nrow(distA.before), inside=nrow(distA.inside), after=nrow(distA.after))
toplotB <- c(before=nrow(distB.before), inside=nrow(distB.inside), after=nrow(distB.after))

toplotAB <- rbind(toplotA, toplotB)
bpAB <- barplot(toplotAB, main="Number of genes that have eqtls at each position relative to itself",
        xlab="Number of Gears", col=c("darkblue","red"),
        legend = c("causal genes", "affected genes"), beside=TRUE, ylim=c(0,max(toplotAB)+50))
text(bpAB, toplotAB, labels=toplotAB, cex=1, pos=3)
