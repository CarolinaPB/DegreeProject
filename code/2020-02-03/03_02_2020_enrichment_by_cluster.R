library(data.table)
data(yeast)


path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))



load(paste0(path,"results/2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))

cluster_table <- data.table(lc$nodeclusters) # nodes(genes) and the cluster they belong to 
clust_w_genes <- list() # list where each element is a list of genes in a cluster (name of list is number of cluster)
for (clusternum in unique(cluster_table$cluster)){
  clust_w_genes[[paste0("cluster_", clusternum)]] <- levels(droplevels(unique(cluster_table[cluster==clusternum]$node)))
}


functional_set = yeast$goslimproc # list containing the gene sets of the GOslim process ontology (Ashburner et al., 2000) for the buddying yeast Saccaromyces Cerevisiae
gene_names <- unique(c(find.effects_TF$geneA, find.effects_TF$geneB)) # all genes in the network

# enrichment of individual clusters with full network
enrichment.res.1 = neat(alist = clust_w_genes, blist = functional_set, network = as.matrix(yeast_net),
            nettype = 'directed', nodes = gene_names, alpha = 0.01, mtc.type = 'fdr')

enrichment.res <- data.table(enrichment.res.1)

enrichment.res.over <- enrichment.res[conclusion=="Overenrichment"] # keep cases where there is overenrichment

plot(lc, type = "graph", node.pies = F, vlabel=T, vshape="circle" , vsize=3, ewidth=1, clusterids = 1)

getNodesIn(lc, clusterids = 1)

cl1.new <- yeast_net[geneA %in% c("YLR264W", "YLR270W", "YLR281C") & geneB=="YGL193C"]


# don't know if it is correct. it's looking for the enrichment of the genes in that cluster in the whole network

# try to check each cluster agains the network of it's own cluster

# network is the subnetwork of this cluster
lc.edges <- data.table(lc$edges)
clust_enrichment <- list()
edges_w_cluster <- merge(lc.edges, yeast2, all=T, by=c("node1", "node2"))
for (num in 1:length(clust_w_genes)){
  gene_names <- clust_w_genes[[num]]
  edges <- edges_w_cluster[cluster==num]
  clust_enrichment[[names(clust_w_genes[num])]] <- data.table(neat(alist = clust_w_genes[num], blist = functional_set, network = as.matrix(edges[,.(node1, node2)]),
                                                    nettype = 'directed', nodes = gene_names, alpha = 0.01, mtc.type = 'fdr'))
}

# No overenriched when using the enrichment by subnetwork
for (el in clust_enrichment){
  if (nrow(el[conclusion=="Overenrichment"])){
    print(el[conclusion=="Overenrichment"])
  }
  
}





#### enrichment of top 10 clusters with higher modularity
load(paste0(respath, "2020-01-29/communityconnectedness.Rdata"))
load(paste0(respath, "2020-01-29/communitymodularity.Rdata"))

con_mod <- data.table(cluster=1:length(cm),connectedness=cconnect, modularity=cm)

con_mod_numnodes <- merge(con_mod,data.table(cluster=as.integer(names(lc$clustsizes)), numnodes=lc$clustsizes))
con_mod_numnodes <- con_mod_numnodes[order(-modularity)]
topmod_clusters <- con_mod_numnodes$cluster[1:10]
topmod_clusters.name <- paste("cluster_",topmod_clusters, sep="")

plot(lc, type = "graph", node.pies = F, vlabel=T, vshape="circle" , vsize=3, ewidth=1, clusterids= topmod_clusters)


enrichment.res.over[A %in% topmod_clusters.name]

getNodesIn(lc, clusterids = 2563)

clust_enrichment[topmod_clusters.name]
