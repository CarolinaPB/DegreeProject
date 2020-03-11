
for (i in 3000:4000){
  gene1 <- coordinates_plot_cor$geneA[i]
  gene2 <- coordinates_plot_cor$geneB[i]
  cor.genes <- rcorr(as.matrix(phenotype[,c(..gene1, ..gene2)]))

  cond <- abs(round(unlist(data.table(cor.genes$r)[1,2]), 6)) == abs(round((unique(coordinates_plot_cor[geneA==gene1 & geneB==gene2])$cor), 6))
  if (!cond){
    print(paste("geneA =", gene1, "geneB =", gene2, "new cor =", round(unlist(data.table(cor.genes$r)[1,2]), 6), "old cor =",round((coordinates_plot_cor[geneA==gene1 & geneB==gene2]$cor), 6)), sep=" ")
  }
}

