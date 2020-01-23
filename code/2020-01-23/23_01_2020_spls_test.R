library(data.table)
library(spls)

path <- "/home/carolpb/DegreeProject/" # use with UPPMAX
respath <- "/proj/snic2019-8-367/private/carol/results/"


phenotype <- fread(paste0(path,"data/SI_Data_01_expressionValues.txt"))
genotype <- fread(paste0(path,"data/SI_Data_03_genotypes.txt"))

cv <- cv.spls(genotype[,2:ncol(genotype)], phenotype[,2:ncol(phenotype)], eta = seq(0.1,0.9,0.1), K = c(5:10) )

save(cv, file=paste0(respath, "2020-01-23/23_01_2020_spls_cv_test.Rdata"))