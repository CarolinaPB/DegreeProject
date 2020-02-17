library("GSEABase")
library("GOstats")
#library(parallel)
library(data.table)
#
# path <- "/Users/Carolina/Documents/GitHub/DegreeProject/"
# respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"

path <- "/home/carolpb/DegreeProject/" # uppmax
respath <- "/proj/snic2019-8-367/private/carol/results/" # uppmax

source(paste0(path, "code/myfunctions.R"))

find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))

# python 10_02_2020_yeastminepy_go_evidencecode.py > /Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-10/yeastmine/yeastmine_evidencecode.txt
ympy.evidence <- fread(paste0(respath, "2020-02-10/yeastmine/yeastmine_evidencecode.txt"))

load(paste0(respath,"2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))
edges <- data.table(lc$edges)

goframeData <- unique(ympy.evidence[,.(GO, evidencecode, gene)])

gs <- getgeneset(goframeData)

#getclusterenrichment(numcluster = 1, geneset = gs, net.cluster = edges)
start_time <- Sys.time()
cl = makeCluster(detectCores(), type="FORK")
clustenrich <- parLapply(cl=cl, unique(edges$cluster), getclusterenrichment, edges, gs)
stopCluster(cl)
clustenrich <- lapply(unique(edges$cluster), getclusterenrichment, edges, gs)
end_time <- Sys.time()
end_time - start_time
save(clustenrich, file=paste0(respath, "2020-02-10/enrichment_bycluster.Rdata"))
print("done")
