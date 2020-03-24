library(data.table)
library(parallel)

path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" # use with own computer
# path <- "/home/carolpb/DegreeProject/" # use with uppmax
#respath <- "/proj/snic2019-8-367/private/carol/results"
# path <- "/home/carolina/DegreeProject/" # use with diprotodon
# respath <- "/home/carolina/DegreeProject/results/"

load(paste0(path, "results/2020-01-10/effectstable.Rdata"))

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


cl = makeCluster(detectCores()-2, type="FORK")
res_table.temp <- parLapply(cl=cl, res, create_res_table)
stopCluster(cl)

res_table <- rbindlist(res_table.temp)
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





# not taking the cor-pval into account
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


#Number of times each gene is the causal one or is on the receiving end
ngeneA <- find.effects_TF[find.effects_TF$`A->B`==T, .N, by=geneA]
ngeneB <- find.effects_TF[find.effects_TF$`B->A`==F, .N, by=geneB]
ngeneA_B <- merge(ngeneA, ngeneB, by.x = "geneA", by.y = "geneB", all=T)
colnames(ngeneA_B) <- c("gene", "A","B")
ngeneA_B[order(-A)]
