library(data.table)
library(parallel)

path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" # use with own computer
# path <- "/home/carolpb/DegreeProject/" # use with uppmax
#respath <- "/proj/snic2019-8-367/private/carol/results"
# path <- "/home/carolina/DegreeProject/" # use with diprotodon
# respath <- "/home/carolina/DegreeProject/results/"

effects_table.cor <- fread(paste0(path, "results/2020-01-10/effectstable.gz"))

source(paste0(path, "code/myfunctions.R"))

nSNPs <- 42052
nGenes <- 5720
snp.pval <- 0.05/(as.numeric(nGenes) * as.numeric(nSNPs))
snp.pval.nsign <- as.numeric(1e-5)
corr.pval <-  0.05/(nGenes*nGenes)



sign_p <- c(1e-17, 1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
non_sign_p <- c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
cor_p <- c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
params0 <- data.table(expand.grid(sign_p=sign_p, non_sign_p=non_sign_p, cor_p=cor_p))
#### Ran on diprotodon

# cl = makeCluster(detectCores()-2, type="FORK")
cl = makeCluster(4, type="FORK")
res <- parApply(cl=cl,params0,1, testparams, effects_table.cor)
stopCluster(cl)

#save(res, file=paste0(respath, "2020-01-20/20_01_2020_resparams.Rdata"))

#####


load(paste0(path,"results/2020-01-20/20_01_2020_resparams.Rdata"))

tab <- data.table(matrix(ncol=9, nrow=length(res)))
names(tab) <- c("res","geneA", "geneB", "eqtl.A", "eqtl.B","unique.genepairs","sign.p", "nonsign.p", "cor.p")

cl = makeCluster(detectCores()-2, type="FORK")
# each row corresponds to one set of parameters
for (n in 1:length(res)){
  tab$res[n] <- n
  tab$geneA[n] <- length(unique(res[[n]]$geneA))
  tab$geneB[n] <- length(unique(res[[n]]$geneB))
  tab$eqtl.A[n] <- length(unique(res[[n]]$eqtl.A))
  tab$eqtl.B[n] <- length(unique(res[[n]]$eqtl.B))
  tab$sign.p[n] <- res[[n]]$sign.p[1]
  tab$nonsign.p[n] <- res[[n]]$non_sign_p[1]
  tab$cor.p[n] <- res[[n]]$cis_p[1]
  tab$unique.genepairs[n] <- res[[n]][(`A->B`==T & `B->A`==F), .N, by=.(geneA, geneB)]
  
  
}
stopCluster(cl)


res_sub <- res[1:2]


create_res_table <- function(tab){
  temp <- data.table(matrix(ncol=8, nrow=1))
  names(temp) <- c("geneA", "geneB", "eqtl.A", "eqtl.B","unique.genepairs","sign.p", "nonsign.p", "cor.p")
  
  temp$geneA <- length(unique(tab$geneA))
  temp$geneB <- length(unique(tab$geneB))
  temp$eqtl.A <- length(unique(tab$eqtl.A))
  temp$eqtl.B <- length(unique(tab$eqtl.B))
  temp$sign.p <- tab$sign.p[1]
  temp$nonsign.p <- tab$non_sign_p[1]
  temp$cor.p <- tab$cis_p[1]
  un.genepairs <- tab[(tab$`A->B`==T & tab$`B->A`==F), .N, by= .(geneA, geneB)]
  temp$unique.genepairs <- nrow(un.genepairs[N==1])
  
  return(temp)
}


cl = makeCluster(detectCores()-2, type="FORK")
res_table.temp <- parLapply(cl=cl, res, create_res_table)
stopCluster(cl)

res_table <- rbindlist(res_table)
res_table <- res_table[order(sign.p, -nonsign.p, cor.p)]

# plot(res_table$sign.p, res_table$unique.genepairs)



par(mfrow=c(2,4))

cor.pvals <- unique(res_table$cor.p)
linetype <- c(1:length(unique(res_table$nonsign.p)))

for (p in cor.pvals){
  xrange <- range(-log(res_table$sign.p)) # set x-axis range
  yrange <- range(res_table$unique.genepairs) # set y-axis range
  plot(xrange, yrange, type = "n", main=paste("corr pval = ", p), xlab = "-log(sign.p)", ylab = " #unique gene pairs") # empty plot
  colors <- rainbow(length(unique(res_table$nonsign.p)))
  for (i in 1:length(unique(res_table$nonsign.p))){
    x <- -log(res_table[nonsign.p==nonsign.p[i] & cor.p==p]$sign.p)
    y <- res_table[nonsign.p==nonsign.p[i] & cor.p==p]$unique.genepairs
    lines(x,y,  type="l", lwd=1.5,lty=linetype[i], col=colors[i])
  }
  #legend(x=20, yrange[2], unique(res_table$nonsign.p), col=colors, lty=linetype, cex=0.8)
}






par(mfrow=c(1,1))

xrange <- range(-log(res_table$sign.p)) # set x-axis range
yrange <- range(res_table$unique.genepairs) # set y-axis range
plot(xrange, yrange, type = "n", xlab = "-log(p)", ylab = "#unique gene pairs",  main="#gene pairs where A->B") # empty plot
colors <- rainbow(length(unique(res_table$nonsign.p)))
linetype <- c(1:length(unique(res_table$nonsign.p)))
for (i in 1:length(unique(res_table$nonsign.p))){
  x <- -log(res_table$sign.p[res_table$nonsign.p==unique(res_table$nonsign.p)[i]])
  y <- res_table$unique.genepairs[res_table$nonsign.p==unique(res_table$nonsign.p)[i]]
  lines(x,y,  type="l", lwd=1.5,lty=linetype[i], col=colors[i])
}
legend(x=20, yrange[2], unique(res_table$nonsign.p), col=colors, lty=linetype, cex=0.8)


