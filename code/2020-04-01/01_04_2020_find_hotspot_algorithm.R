#rules
# if runs of F on both sides longer than something
# and the fraction of genes in between with count.A > x is greater than something
# then call it a hotspot

causalgenes.pos.count1 <- unique(coordinates_plot_links[order(chr.A,start.A)][,.(geneA, start.A, end.A, chr.A, count.A)])

allgenes.location <- genepos

## transform chromosome ids into numbers
# remove "chr" part of the chromosome name
allgenes.location$chr.id <- gsub('chr', '', allgenes.location$chr.id)
# convert roman chromosome numbers to numbers
allgenes.location$chr.id <- as.numeric(as.roman(allgenes.location$chr.id))

causalgenes.pos.count <- merge(causalgenes.pos.count1, allgenes.location, by.x=c("geneA", "chr.A"), by.y=c("gene", "chr.id"))
causalgenes.pos.count <- causalgenes.pos.count[order(chr.A, start.A, end.A)]

rle.res.list <- list()
for (chr in unique(causalgenes.pos.count$chr.A)){
  cond <- causalgenes.pos.count[chr.A == chr]$count.A >= 10
  rle.res <- rle(cond)
  rle.res.dt <- data.table(idx=1:length(rle.res$lengths),lengths=rle.res$lengths, values=rle.res$values)
  rle.res.list[[chr]] <- rle.res.dt
}

# rle.res.list - list of tables with results from rle. one table per chromosome
# True values represent pre-hubs or hotspots
for (el in 1:length(rle.res.list)){ # iterate through each chromosome's table
 tab <- rle.res.list[[el]]
 for (i in tab[values==T]$idx){ # iterate through all potential hubs/hotspots
    if (i == 1){
      # hotspot <- causalgenes.pos.count[chr.A==el][1:tab[idx==i]$lengths,]
    } else {
     previous <- i-1
     after <- i+1
     if (tab[idx==previous]$lengths > 2 & tab[idx==after]$lengths > 2){
       print(paste(el, i,"is a hotspot"))
     }
     
   }
 }
}

chr <- 7
id <- 2
id <- 12
id <- 4

tab <- rle.res.list[[chr]]
id.group <-vector()
id.keep <- vector()
for (id in tab[values==T & lengths>=3]$idx){
  if (tab[idx==id-1]$lengths > 2 & tab[idx==id+1]$lengths > 2){
    id.keep <- c(id.keep, id)
  } else if (tab[idx==id+1]$lengths < 3){
    id.group <- c(id.group, id, id+1)
  }
}

id <- 7
getrows <- (sum(tab[idx<=id-1]$lengths)+1):sum(tab[idx<=id]$lengths)
causalgenes.pos.count[chr.A==chr][getrows,]

plot_sorted_coordinates(coordinates_plot_cor[chr.A==12], separator = separator, col = coordinates_plot_cor[chr.A==chr]$col)

tab <- rle.res.list[[chr]]
test <- tab[(lengths<3 & values==F) | (values==T)]


res <- list()
for (el in 1:length(rle.res.list)){ # iterate through each chromosome's table
  tab <- rle.res.list[[el]]
  keep <- vector()
  if (nrow(tab) >1){
    for (numrow in 1:(nrow(tab) - 1)) {
      if (tab[numrow, ]$values == T & tab[numrow, ]$lengths >3) {
        keep <- c(keep, tab[numrow, ]$idx)
      }
    }
  } else if(nrow(tab) ==1 & tab$values==T){
    keep <- c(keep, tab[numrow, ]$idx)
  }
  res[[el]] <- keep 
}

# for each consecutive pair, check if the length for idx between them is < 3

tt <- res[[7]]
tab <- rle.res.list[[7]]
group <- vector()
for (i in 1:(length(tt)-1)){
  if (tt[i+1]-tt[i] ==2){
    if(nrow(tab[idx==(tt[i]+1) & lengths <3])>0){
      group <- c(tt[i], tt[i+1])
    }
  }
}

