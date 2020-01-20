library(data.table)
library(parallel)

#path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" # use with own computer
# path <- "/home/carolpb/DegreeProject/" # use with uppmax
#respath <- "/proj/snic2019-8-367/private/carol/results"
path <- "/home/carolina/DegreeProject/" # use with diprotodon
respath <- "/home/carolina/DegreeProject/results/"

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


# cl = makeCluster(detectCores()-2, type="FORK")
cl = makeCluster(4, type="FORK")
res <- parApply(cl=cl,params0,1, testparams, effects_table.cor)
stopCluster(cl)

save(res, file=paste0(respath, "2020-01-20/20_01_2020_resparams.Rdata"))