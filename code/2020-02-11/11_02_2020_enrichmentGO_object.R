library("GSEABase")
library("GOstats")
library(parallel)
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


getgeneset <- function(df){
  # create gene set to use for enrichment
  # requires a dataframe with three columns
  # first column is GO ids
  # second column is evidence codes
  # third column is gene ids that match those GO terms
  
  goFrame=GOFrame(goframeData,organism="Saccharomyces cerevisiae")
  goAllFrame=GOAllFrame(goFrame)
  
  return(GeneSetCollection(goAllFrame, setType = GOCollection()))
}

getclusterenrichment <- function(numcluster, net.cluster, geneset){
  # get enrichment for each cluster
  # requires a table with three columns
  # first column is the causal genes
  # second column is the affected genes
  # third column is the cluster number they belong to
  # numcluster - number of the cluster to analyse
  # geneset - a geneset created from the GO terms, genes and evidence code - result from getgeneset
  
  genes <- levels(droplevels(net.cluster[cluster==numcluster]$node1))
  universe = unique(c(levels(droplevels(net.cluster$node1)), levels(droplevels(net.cluster$node2))))
  params <- GSEAGOHyperGParams(name="first try",
                               geneSetCollection=geneset,
                               geneIds = genes,
                               universeGeneIds = universe,
                               ontology = "BP",
                               pvalueCutoff = 0.05,
                               conditional = FALSE,
                               testDirection = "over")
  Over <- hyperGTest(params)
  # res_over <- data.table(summary(Over))
  # return(res_over)
  return(Over)
}


load(paste0(respath,"2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))
edges <- data.table(lc$edges)

goframeData <- unique(ympy.evidence[,.(GO, evidencecode, gene)])

gs <- getgeneset(goframeData)

#getclusterenrichment(numcluster = 1, geneset = gs, net.cluster = edges)
start_time <- Sys.time()
cl = makeCluster(detectCores(), type="FORK")
clustenrich <- parLapply(cl=cl, unique(edges$cluster), getclusterenrichment, edges, gs)
stopCluster(cl)
end_time <- Sys.time()
end_time - start_time
save(clustenrich, file=paste0(respath, "2020-02-11/enrichment_bycluster_object.Rdata"))
print("done")

