library(data.table)
library(linkcomm)

path <- "/home/carolpb/DegreeProject/" # use with uppmax
respath <- "/proj/snic2019-8-367/private/carol/results/"

find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))
load(paste0(respath,"2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))


pdf(paste0(respath, "2020-01-28/lc.summary.pdf"))
plot(lc, type = "summary")
dev.off()

