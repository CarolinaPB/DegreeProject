library(data.table)
data(yeast)


path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))



load(paste0(path,"results/2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))

cluster_table <- data.table(lc$nodeclusters)
clust_w_genes <- list()
for (clusternum in unique(cluster_table$cluster)){
  clust_w_genes[[paste0("cluster_", clusternum)]] <- levels(droplevels(unique(cluster_table[cluster==clusternum]$node)))
}


functional_set = yeast$goslimproc
yeast_net <- find.effects_TF[,.(geneA, geneB)]
gene_names <- unique(c(find.effects_TF$geneA, find.effects_TF$geneB))


enrichment.res.1 = neat(alist = clust_w_genes, blist = functional_set, network = as.matrix(yeast_net),
            nettype = 'directed', nodes = gene_names, alpha = 0.01, mtc.type = 'fdr')

enrichment.res <- data.table(enrichment.res.1)

enrichment.res.over <- enrichment.res[conclusion=="Overenrichment"]

plot(lc, type = "graph", node.pies = F, vlabel=T, vshape="circle" , vsize=3, ewidth=1, clusterids = 1)

getNodesIn(lc, clusterids = 1)

cl1.new <- yeast_net[geneA %in% c("YLR264W", "YLR270W", "YLR281C") & geneB=="YGL193C"]


# don't know if it is correct. it's looking for the enrichment of the genes in that cluster in the whole network

# try to check each cluster agains the network of it's own cluster

