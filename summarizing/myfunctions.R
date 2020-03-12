# Functions I've created for my analysis + some useful ones found on the internet


effect_eqtl_gene <- function(res, pheno=phenotype, geno=genotype){
  # function that calculates the effect of a qtl on a gene - using anova
  # outputs a string with anova.pval__anova.r2
  # res - data.table with 2 columns;
  # first column: gene
  # second column: eqtl
  lmp <- function (modelobject) {
    # function to get the p-value out of a lm
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
  }

  gene <- res[1]
  eqtl <- res[2]

  anv <- lm(unlist(pheno[ , ..gene]) ~ unlist(geno[ ,..eqtl]))
  anv.res <- paste(lmp(anv), summary(anv)$adj.r.squared, sep="__")

  return(anv.res)
}

merge_after_anova <- function(res.effect_eqtl_gene, gene.AB, eqtl.AB, effects.table){
  # merge the results of the anova with the main table
  # requires library(tidyr) for the separate
  gene.AB <- toupper(gene.AB)
  eqtl.AB <- toupper(eqtl.AB)

  pval.name <- paste0("eqtl",eqtl.AB, "_", "gene", gene.AB, ".pval")
  r2.name <- paste0("eqtl",eqtl.AB, "_", "gene", gene.AB, ".r2")

  res.sep <- res.effect_eqtl_gene %>% separate(anv.res, c(pval.name, r2.name), sep="__")

  setnames(res.sep,"gene", paste0("gene", gene.AB))
  setnames(res.sep,"eqtl", paste0("eqtl.", eqtl.AB))

  effects.eqtlB <- merge(effects.table, res.sep, by=c(paste0("gene", gene.AB),paste0("eqtl.", eqtl.AB)), all.x=T)

  return(effects.eqtlB)
}

get_num_genos <- function(res, genotype){
  # get the number of each geno per snp
  # res is a data.table with one row - the snp ids
  snp <- res[1]

  geno_1 <- sum(ifelse(genotype[,..snp] == 1, 1,0))
  geno_neg1 <- sum(ifelse(genotype[,..snp] == -1, 1,0))

  return(paste(geno_1, geno_neg1, sep="__"))
}

flat_cor_mat <- function(cor_r, cor_p){
  #This function provides a simple formatting of a correlation matrix
  #into a table with 4 columns containing :
  # Column 1 : row names (variable 1 for the correlation test)
  # Column 2 : column names (variable 2 for the correlation test)
  # Column 3 : the correlation coefficients
  # Column 4 : the p-values of the correlations

  library("dplyr")
  library(tidyr)
  library(tibble)
  cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
  cor_r <- gather(cor_r, column, cor, -1)
  cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
  cor_p <- gather(cor_p, column, p, -1)
  cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
}

expand.grid.faster <- function(seq1,seq2) {
  # faster alternative to expand.grid
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}

create_ini_table <- function(eqtl.table, genesB, var.exp.lim){
  # requires data.table, tidyr
  # takes a table with 4 columns: gene, eqtl, cis and var.exp and creates a table with all the possible combinations
  # of geneA-eqtlA with geneB-eqtlB
  # cis is a True/false - if the eqtl is in cis with the gene
  # genesB is the genes we want to compare the gene-eqtl pairs with
  # var.exp.lim is chosen by the user
  # For the first set of pairs, only the gene-eqtl that are in cis and where the variance explained is
  # above the set limit will be kept

  # Keep only the gene-eqtl pairs where the eqtl is in cis with the gene and where the var.exp >0.1
  eqtl.tableA <- eqtl.table[cis==T & var.exp > var.exp.lim]

  # unite columns so they act as one block of information
  eqtl.tableA.unite <- eqtl.tableA %>% unite(infoA, gene, pmarker, cis, var.exp, sep = "__")

  # get table with all genes that are going to be tested with eqtlA (geneB migth have an eqtl or not)
  eqtl.tableB<- merge(data.table(genesB),eqtl.table, by.x="genesB", by.y="gene", all.x=T)

  # unite columns so they act as one block of information
  eqtl.tableB.unite <- eqtl.tableB %>% unite(infoB, genesB, pmarker, cis, var.exp, sep = "__")

  # get all combinations of geneA and geneB with corresponding eqtls
  eqtl.table.combineAB <- expand.grid.faster(eqtl.tableA.unite$infoA,eqtl.tableB.unite$infoB)
  eqtl.table.combineAB.dt <- data.table(eqtl.table.combineAB)
  setnames(eqtl.table.combineAB.dt, old=c("Var1", "Var2"), new=c("infoA", "infoB"))

  # separate info blocks into normal columns again
  eqtl.table.sepA <- eqtl.table.combineAB.dt %>% separate(infoA, c("geneA", "eqtl.A", "cis.A", "var.exp.A"), sep="__")
  eqtl.table.sepAB <- eqtl.table.sepA %>% separate(infoB, c("geneB", "eqtl.B", "cis.B", "var.exp.B"), sep="__")
}

find.effects_fun <- function(effects.table, snp.pval, snp.pval.nsign){
  # the input table must have two columns called "eqtlA_geneB.pval" and "eqtlB_geneA.pval", which contain the p-val of the anova of eqtlX on geneY
  # This function does not take into account if the eqtl is in cis with the gene, so the input table must have been subseted to only contain eqtl-gene pairs that are in cis

  # adds two columns to the input table: A->B and B->A.
  # These columns can have a value of
  #   * NA- if with the chosen p-val cutoffs we can't say if X affects Y
  #   * TRUE - if X affects Y
  #   * FALSE - if X does not affect Y

  # if both genes have same eqtl, then it's the third case, C, where both are affected by the same.
  # we can't say anything about the relationship between A and B.
  # We just know that a third factor is affecting them

  effects.table[, ("A->B") := NA]; effects.table[as.numeric(eqtlA_geneB.pval) < snp.pval, ("A->B") := T];
  #effects.table[(as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) ), ("A->B") := F]
  effects.table[(as.numeric(eqtlA_geneB.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B))), ("A->B") := F]

  effects.table[, ("B->A") := NA]; effects.table[as.numeric(eqtlB_geneA.pval) < snp.pval, ("B->A") := T];
  #effects.table[(as.numeric(eqtlB_geneA.pval) > as.numeric(snp.pval.nsign)), ("B->A") := F]
  effects.table[(as.numeric(eqtlB_geneA.pval) > as.numeric(snp.pval.nsign) & (as.character(eqtl.A) != as.character(eqtl.B))), ("B->A") := F]

  return(effects.table)
}

testparams <- function(params, tab, var.exp){
  #function to get cases where T and F with different params
  # to run in parallel
  sign_p <- params[1]
  non_sign_p <- params[2]
  cor_p <- params[3]


  temp <- tab[cor.pval < cor_p & cis.A ==T & cis.B==T & geneA!=geneB]
  temp <- find.effects_fun(temp, sign_p, non_sign_p)

  temp_TF.1 <- temp[temp$`A->B`==T & temp$`B->A`==F & var.exp.A > var.exp, .(geneA, geneB, eqtl.A, eqtl.B, `A->B`, `B->A`)]
  temp_TF.2 <- rbind(temp_TF.1, temp[temp$`A->B`==F & temp$`B->A`==T & var.exp.B > var.exp, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
  temp_TF <- unique(temp_TF.2)

  temp_TF$sign.p <- sign_p
  temp_TF$non_sign.p <- non_sign_p
  temp_TF$cor.p <- cor_p

  return(temp_TF)
}

create_res_table <- function(tab){
  # create table with counts for several categories (see below)
  # will need to be "collapsed with rbindlist()
  
  temp <- data.table(matrix(ncol=10, nrow=1))
  names(temp) <- c("geneA", "geneB", "eqtl.A", "eqtl.B","unique.genepairs", "unique.geneeqtlpairs.A","unique.geneeqtlpairs.B","sign.p", "nonsign.p", "cor.p")

  temp$geneA <- length(unique(tab$geneA))   # number of unique genesA
  temp$geneB <- length(unique(tab$geneB))   # number of unique genesB
  temp$eqtl.A <- length(unique(tab$eqtl.A)) # number of unique eqtlsA
  temp$eqtl.B <- length(unique(tab$eqtl.B)) # number of unique eqtlsB
  temp$sign.p <- tab$sign.p[1]              # sign.pval cutoff
  temp$nonsign.p <- tab$non_sign_p[1]       # non sign.pval cutoff
  temp$cor.p <- tab$cis_p[1]                # corr.pval cutoff
  un.genepairs <- tab[(tab$`A->B`==T & tab$`B->A`==F), .N, by= .(geneA, geneB)] 
  temp$unique.genepairs <- nrow(un.genepairs[N==1])   # number of unique gene pairs
  un.geneeqtl.A <- tab[(tab$`A->B`==T & tab$`B->A`==F), .N, by= .(geneA, eqtl.A)]   
  un.geneeqtl.B <- tab[(tab$`A->B`==T & tab$`B->A`==F), .N, by= .(geneB, eqtl.B)]   
  temp$unique.geneeqtlpairs.A <- nrow(un.geneeqtl.A)      # number of geneA-eqtlA pairs
  temp$unique.geneeqtlpairs.B <- nrow(un.geneeqtl.B)      # number of geneB-eqtlB pairs

  return(temp)
}

getgeneset <- function(goframeData){
  # library("GSEABase")
  # library("GOstats")

  # create gene set to use for enrichment
  # requires a dataframe with three columns
  # first column is GO ids
  # second column is evidence codes
  # third column is gene ids that match those GO terms

  goFrame=GOFrame(goframeData,organism="Saccharomyces cerevisiae")
  goAllFrame=GOAllFrame(goFrame)

  return(GeneSetCollection(goAllFrame, setType = GOCollection()))
}

getclusterenrichment <- function(numcluster, net.cluster, geneset){
  # library("GSEABase")
  # library("GOstats")

  # get enrichment for each cluster
  # requires a table with three columns
  # first column is the causal genes
  # second column is the affected genes
  # third column is the cluster number they belong to
  # numcluster - number of the cluster to analyse
  # geneset - a geneset created from the GO terms, genes and evidence code - result from getgeneset

  genes <- levels(droplevels(net.cluster[cluster==numcluster]$node1))
  universe = unique(c(levels(droplevels(net.cluster$node1)), levels(droplevels(net.cluster$node2))))
  params <- GSEAGOHyperGParams(name="first try",
                               geneSetCollection=geneset,
                               geneIds = genes,
                               universeGeneIds = universe,
                               ontology = "BP",
                               pvalueCutoff = 0.05,
                               conditional = FALSE,
                               testDirection = "over")
  Over <- hyperGTest(params)
  res_over <- data.table(summary(Over))
  return(res_over)
}

getenrichment <- function(geneset, universe, interestinggenes, pval=0.05, cond = FALSE){
  # need these libraries

  # library("GSEABase")
  # library("GOstats")

  # get enrichment for a set of genes
  # geneset - a geneset created from the GO terms, genes and evidence code - result from getgeneset
  # geneAB can be "geneA" or "geneB"
  # returns an object with the results of the enrichment analysis.
  # summary table can be created with summary(results)
  # info can be extracted with the functions found here ?geneCounts

  # genes <- unlist(unique(find.effects_TF[,..geneAB]))
  # universe <- unique(c(find.effects_TF$geneA, find.effects_TF$geneB))

  params <- GSEAGOHyperGParams(name="GO enrichment",
                               geneSetCollection=geneset,
                               geneIds = interestinggenes,
                               universeGeneIds = universe,
                               ontology = "BP",
                               pvalueCutoff = pval,
                               conditional = cond,
                               testDirection = "over")
  Over <- hyperGTest(params)
}

plot_enrichment_heatmap <- function(dt, ...){
  # accepts data.table where the first column is the GO term, the second is the causal genes and the third is the affected genes
  # plots a heatmap with custom color in a log10 scale
  
  paletteLength <- 50
  myColor <- colorRampPalette(c("white", "darkblue","blue","lightblue","green","lightyellow", "yellow","orange", "red", "darkred"))(paletteLength)
  setnames(dt, c("term", "Causal", "Affected"))
  myBreaks <- c(seq(min(dt[,c(2,3)]), max(dt[,c(2, 3)])/paletteLength-0.000001, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(dt[,c(2,3)])/paletteLength, max(dt[,c(2, 3)]), length.out=floor(paletteLength/2)))
  
  pheatmap(as.matrix(dt, rownames = 1), 
           fontsize_row=3, cellwidth=70, cluster_col=F, cluster_rows=T, 
           fontsize = 5, treeheight_row = 0, color=myColor, breaks=myBreaks, border_color="white", 
           main="Heatmap enrichment p-value (-log10)", ...)
  
}

sort_by_chr <- function(vchr, genepairs_pos, separator){
  # Sorts genes' chromosome positions to be plotted
  # takes a numeric vector of chromosome numbers (should be ordered by the order you want to plot by)
  # takes a table where there must be the following columns:
  # geneA  - to be plotted on the x axis
  # geneB  - to be plotted on the y axis
  # start.A - position of geneA
  # start.B - position of geneB
  # chr.A - chromosome of geneA
  # chr.B - chromosome of geneB
  # interval between genes (for plotting purposes)
  
  res <- data.table(genepairs_pos)
  
  for (i in 1:length(vchr)){
    if (vchr[i] != 1) {
      previous <- vchr[i-1]
      res[chr.A == vchr[i]]$start.A <- max(res[chr.A==previous]$start.A) + res[chr.A == vchr[i]]$start.A + separator
      res[chr.B == vchr[i]]$start.B <- max(res[chr.B==previous]$start.B) + res[chr.B == vchr[i]]$start.B + separator
    }
  }
  return(res)
}


plot_sorted_coordinates <- function(coordinates_plot, separator, ...){
  # plots the gene position coordinates by chromosome using the result table from sort_by_chr()
  # accepts extra parameters for the plot function
  
  # takes a table where there must be the following columns:
  # geneA  - to be plotted on the x axis
  # geneB  - to be plotted on the y axis
  # start.A - position of geneA
  # start.B - position of geneB
  # chr.A - chromosome of geneA
  # chr.B - chromosome of geneB
  # interval between genes (for plotting purposes)
  
  # the separator should be the same as the one used for the sort_by_chr() function
  
  plot(coordinates_plot$start.A, coordinates_plot$start.B, pch=".", axes=F, 
       xlab = "Causal gene position (chr)", ylab = "Affected gene position (chr)", ...)
  
  chrA <- unique(coordinates_plot[order(chr.A)]$chr.A)
  chrB <- unique(coordinates_plot[order(chr.B)]$chr.B)
  # nchr <- max(coordinates_plot$chr.A)
  # add chromosome separators
  for (ch in chrA){
    abline(v= max(coordinates_plot[chr.A==ch]$start.A)+separator/2, col="lightblue", lty=2) 
  }
  for (ch in chrB){
    abline(h= max(coordinates_plot[chr.B==ch]$start.B)+separator/2, col="lightblue", lty=2) 
  }
  # add x and y axis chromosomes to the axis
  axis(1, at=sapply(chrA, function(i){min(coordinates_plot[chr.A==i]$start.A) + (max(coordinates_plot[chr.A==i]$start.A) - min(coordinates_plot[chr.A==i]$start.A))/2}), labels=as.roman(chrA), tick=FALSE)
  axis(2, at=sapply(chrB, function(i){min(coordinates_plot[chr.B==i]$start.B) + (max(coordinates_plot[chr.B==i]$start.B) - min(coordinates_plot[chr.B==i]$start.B))/2}), labels=as.roman(chrB), tick=FALSE)
}


genes_inside_hotspot <- function(hotspot_data, positions_table, chromosome, leftlim=3, rightlim=4){
  # function that takes a table with chromosome number, and hotspot left and 
  # right intervals and finds the genes in my dataset that are inside those intervals
  # prints the genes in the chosen chromosomes that are in the hotspots
  
  # hotspot_data: must have at least three columns: 
  #   chr - chromosome number
  #   bootstrapIntervalLeft - left limit position
  #   bootstrapIntervalRigth - right limit position
  # positions_table - must have at least three columns:
  #   gene names in the first column
  #   chr.A - chromosome number
  #   start.A - gene start position
  # chromosome - chromosome vector of the chromosomes whose genes should be tested
  # lefttlim/righttlim - number of column to be used at the left and right limit of the hotspot
  
  total_count <- 0
  res <- list()
  indx <- 1
  for (chr in chromosome){
    count <- 0
    print(paste("genes in chromosome", chr, "that are in the described hotspots", sep=" "))
    for (num in 1:nrow(hotspot_data[chromosome==chr])){
      left <- hotspot_data[chromosome==chr][num,..leftlim]
      right <- hotspot_data[chromosome==chr][num, ..rightlim]
      
      gene_inside_hotspot <- unlist(unique(positions_table[chr.A==chr][between(start.A, unlist(left), unlist(right))][,1]))
      if (length(gene_inside_hotspot) >0){
        print(unname(gene_inside_hotspot))
        count <- count+1
        total_count <- total_count +length(gene_inside_hotspot)
        res[[indx]] <- unname(gene_inside_hotspot)
        names(res)[indx] <- hotspot_data[chromosome==chr][num,]$hotspotMarker
        indx <- indx+1
      }
    }
    print(paste("your genes overlap with",count, "hotspots", sep=" "))
    cat("\n")
  }
  print(paste("In total there are", total_count, "genes that overlap with the hotspots", sep=" "))
  return(res)
}

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