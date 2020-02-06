library(data.table)
library(linkcomm)
path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
source(paste0(path, "code/myfunctions.R"))


#   ____________________________________________________________________________
#   Parameters                                                              ####

# var.exp.lim <- 0.1
# nSNPs <- 42052
# nGenes <- 5720
# snp.pval <- 0.01
# snp.pval.nsign <- as.numeric(1e-5)
# corr.pval <- 0.05/choose(nGenes,2)


#   ____________________________________________________________________________
#   Load files                                                              ####
find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))
load(paste0(path,"results/2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))

#fwrite(find.effects_TF[,.(geneA, geneB)], paste0(respath, "2020-01-27/edgelist.csv"), na=NA)
#   ____________________________________________________________________________
#   Tests                                                                   ####

print(lc)

# A data frame consisting of 2 columns; 
# the first contains node names, and the second contains single community IDs for each node. 
# head(lc$clusters)All communities and their nodes are represented, but not necessarily all nodes
head(lc$nodeclusters)

# A list of integer vectors containing the link IDs that belong to each community. 
# each list has the members of that community
# Community IDs are the numerical position of the communities in the list.
head(lc$clusters)


# A data frame with 3 columns; 
# the first two contain genes(nodes) that interact with each other
# the third is an integer vector of community IDs indicating community membership for each link.

lc.edges <- data.table(lc$edges)
lc.edges


# Number of communities each gene(node) belongs to

head(lc$numclusters)


# Names are community IDs
# integer values indicate the number of nodes that belong in each community.

head(lc$clustsizes)
tail(lc$clustsizes)

# Nodes that belong to a lot of different communities will get the largest community centrality scores
# nodes that belong to overlapping, nested or few communities will get the lowest scores.
#The 5 most central genes
# cc <- getCommunityCentrality(lc)
#save(cc, file=paste0(respath, "2020-01-29/communitycentrality.Rdata"))
load(paste0(respath, "2020-01-29/communitycentrality.Rdata"))

head(sort(cc, decreasing = TRUE))

# plot genes that belong to the most communities
# pdf(lc, file = paste0(respath, "2020-01-29/communitymembership.pdf"))
plot(lc, type = "members")
# dev.off()
# Numbers on the right - number of communities a gene belongs to
# Numbers on the top - community(cluster) number

# plot a specific cluster (or several )
gr <- graph_from_edgelist(as.matrix(lc.edges[cluster %in% c(974, 988, 964) , .(node1, node2)]), directed=T)
plot(gr, layout=layout_with_graphopt, edge.arrow.size=.4, vertex.label=NA, edge.curved=.1)

gr <- graph_from_edgelist(as.matrix(lc.edges[cluster==974 , .(node1, node2)]), directed=T)
plot(gr, layout=layout_with_graphopt, edge.arrow.size=.4, vertex.label=NA, edge.curved=.1)

# two communities that share genes
gr <- graph_from_edgelist(as.matrix(lc.edges[cluster %in% c(1, 7) , .(node1, node2)]), directed=T)
plot(gr, layout=layout_with_graphopt, edge.arrow.size=.4, edge.curved=.1)
plot(lc, type = "graph", node.pies = T, vlabel=F, clusterids = c(1, 7))

# get shared genes between two communities

get.shared.nodes(lc, comms = c(1,7))

# get nodes in certain communities
getNodesIn(lc, clusterids = c(1,7))


# returns the meta-communities (clusters of community IDs) (when you cut the dendogram at different heights)

cr.20 <- getClusterRelatedness(lc, hcmethod = "ward", cutat = 20)
cr.5 <- getClusterRelatedness(lc, hcmethod = "ward", cutat = 5)
getClusterRelatedness(lc, hcmethod = "ward")


cr.pdmax <- getClusterRelatedness(lc, hcmethod = "ward", cutat = lc$pdmax)


# Networks with high modularity have dense connections between the nodes within modules 
# but sparse connections between nodes in different modules. 
# Modularity is often used in optimization methods for detecting community structure in networks. 

# measure of how relatively outwardly or inwardly connected a community is.
# names are community IDs and the numbers are community connectedness or modularity scores.
#cm <- getCommunityConnectedness(lc, conn = "modularity")
#save(cm, file=paste0(respath, "2020-01-29/communityconnectedness.Rdata"))
load(paste0(respath, "2020-01-29/communityconnectedness.Rdata"))
plot(lc, type = "commsumm", summary = "modularity")



# 10 communities with highest modularity - dense connections within community and sparce between nodes in different communities
# pdf(paste0(respath, "2020-01-29/top10modularity_barplot.pdf"))
plot(lc, type = "commsumm", summary = "modularity", clusterids=as.integer(names(sort(cm, decreasing = T)[1:10])))
# dev.off()
# pdf(paste0(respath, "2020-01-29/top10modularity_network.pdf"))
plot(lc, type = "graph", node.pies = T, vlabel=T, clusterids = as.integer(names(sort(cm, decreasing = T)[1:10])))
# dev.off()


