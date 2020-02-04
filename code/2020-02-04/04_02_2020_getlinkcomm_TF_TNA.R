library(data.table)
library(linkcomm)

#path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
path <- "/home/carolpb/DegreeProject/" # use with uppmax
respath <- "/proj/snic2019-8-367/private/carol/results"

find.effects_TF_TNA <- fread(paste0(respath, "2020-02-04/findeffects_TF_TNA.gz"))

lc.TF_TNA <- getLinkCommunities(find.effects_TF_TNA[,.(geneA, geneB)], directed = T)


save(lc.TF_TNA, file=paste0(respath,"2020-02-04/04_02_2020_linkedcommunities_TF_TNA.Rdata"))