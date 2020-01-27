library(data.table)
library(linkcomm)

#path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
path <- "/home/carolpb/DegreeProject/" # use with uppmax
respath <- "/proj/snic2019-8-367/private/carol/results"

source(paste0(path, "code/myfunctions.R"))

##### PARAMETERS #####
var.exp.lim <- 0.1

# nSNPs <- length(colnames(genotype))-1
# nGenes <- length(colnames(phenotype))-1

nSNPs <- 42052
nGenes <- 5720

snp.pval <- 0.01
snp.pval.nsign <- as.numeric(1e-5)

corr.pval <- 0.05/choose(nGenes,2)


######
effects_table.cor <- fread(paste0(path, "results/2020-01-10/effectstable.gz"))

find.effects <- effects_table.cor[cor.pval < corr.pval & cis.A ==T & cis.B==T & geneA!=geneB]

find.effects <- find.effects_fun(find.effects, snp.pval, snp.pval.nsign)

####
find.effects_TF.1 <- find.effects[find.effects$`A->B`==T & find.effects$`B->A`==F, .(geneA, geneB, eqtl.A, eqtl.B, `A->B`, `B->A`)]
find.effects_TF.2 <- rbind(find.effects_TF.1, find.effects[find.effects$`A->B`==F & find.effects$`B->A`==T, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
find.effects_TF <- unique(find.effects_TF.2)

#### communities stuff
lc <- getLinkCommunities(find.effects_TF[,.(geneA, geneB)], directed = T)

save(lc, file=paste0(respath,"2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))