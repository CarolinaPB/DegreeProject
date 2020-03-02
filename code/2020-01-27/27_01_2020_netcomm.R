library(data.table)
library(linkcomm)

#path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
path <- "/home/carolpb/DegreeProject/" # use with uppmax
respath <- "/proj/snic2019-8-367/private/carol/results"

source(paste0(path, "code/myfunctions.R"))
find.effects_TF <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-03-02/findeffects_TF_newparams.gz")

#### communities stuff
lc <- getLinkCommunities(find.effects_TF[,.(geneA, geneB)], directed = T)

save(lc, file=paste0(respath,"2020-03-02/02_03_2020_linkedcommunities_TF.Rdata"))
