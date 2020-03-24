library(data.table)
path <- "/Users/Carolina/Documents/GitHub/DegreeProject/"
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
load(paste0(path, "results/2020-01-10/effectstable.Rdata"))

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
find.effects <- effects_table.cor[cor.pval < corr.pval & cis.A ==T & cis.B==T & geneA!=geneB]

find.effects <- find.effects_fun(find.effects, snp.pval, snp.pval.nsign)

####
find.effects_TF_TNA.1 <- find.effects[(find.effects$`A->B`==T & find.effects$`B->A`==F), .(geneA, geneB, eqtl.A, eqtl.B, `A->B`, `B->A`)]
find.effects_TF_TNA.2 <- rbind(find.effects_TF_TNA.1, find.effects[find.effects$`A->B`==F & find.effects$`B->A`==T, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
find.effects_TF_TNA.3 <- rbind(find.effects_TF_TNA.2, find.effects[find.effects$`A->B`==T & is.na(find.effects$`B->A`), .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])

find.effects_TF_TNA.4 <- rbind(find.effects_TF_TNA.3, find.effects[is.na(find.effects$`A->B`) & find.effects$`B->A`==T, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
find.effects_TF_TNA <- unique(find.effects_TF_TNA.4)


fwrite(find.effects_TF_TNA, paste0(respath, "2020-02-04/findeffects_TF_TNA.gz"), na = NA)
find.effects_TF_TNA <- fread(paste0(respath, "2020-02-04/findeffects_TF_TNA.gz"))
find.effects_TF_TNA[geneB=="YGL193C" & (geneA %in% c("YLR264W", "YLR270W", "YLR281C"))]




lc.TF_TNA <- getLinkCommunities(find.effects_TF_TNA[,.(geneA, geneB)], directed = T)
