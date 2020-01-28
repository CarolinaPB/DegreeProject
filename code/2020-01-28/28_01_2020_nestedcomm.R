library(data.table)
library(linkcomm)

path <- "/home/carolpb/DegreeProject/" # use with uppmax
respath <- "/proj/snic2019-8-367/private/carol/results"

find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))
load(paste0(respath,"2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))


nested <- getAllNestedComm(lc)

save(nested, file=paste0(respath, "2020-01-28/28_01_2020_nestedcomm.Rdata"))