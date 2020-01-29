library(data.table)
library(linkcomm)

path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"

source("code/myfunctions.R")

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
effects_table.cor <- fread(paste0(path, "results/2020-01-10/effectstable.gz"))

find.effects <- effects_table.cor[cor.pval < corr.pval & cis.A ==T & cis.B==T & geneA!=geneB]

find.effects <- find.effects_fun(find.effects, snp.pval, snp.pval.nsign)

####
find.effects_TF.1 <- find.effects[find.effects$`A->B`==T & find.effects$`B->A`==F, .(geneA, geneB, eqtl.A, eqtl.B, `A->B`, `B->A`)]
find.effects_TF.2 <- rbind(find.effects_TF.1, find.effects[find.effects$`A->B`==F & find.effects$`B->A`==T, .(geneA=geneB, geneB=geneA, eqtl.A=eqtl.B, eqtl.B=eqtl.A, `A->B`=`B->A`, `B->A`=`A->B`)])
find.effects_TF <- unique(find.effects_TF.2)

fwrite(find.effects_TF, paste0(respath, "2020-01-27/findeffects_TF.gz"))
find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))

#### communities stuff
# lc <- getLinkCommunities(find.effects_TF[,.(geneA, geneB)], directed = T)
# 
# save(lc, file=paste0(path,"results/2020-01-24/24_01_2020_linkedcommunities_TF.Rdata"))

load(paste0(path,"results/2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))

print(lc)


lc$nodeclusters
lc$edges
lc$numclusters
lc$clustsizes

cr <- getClusterRelatedness(lc, hcmethod = "ward")
cutdend <- cutDendrogramAt(cr,lc=lc, cutat = 20)

lc.edges <- data.table(lc$edges)
gr <- graph_from_edgelist(as.matrix(lc.edges[cluster==3040, .(node1, node2)]), directed=F)

plot(gr, layout=layout_with_graphopt)




plot(lc, type = "graph", layout = layout.fruchterman.reingold, vlabel = FALSE, node.pies=F)

# The following displays only the nodes that belong to 3 or more communities:
#plot(lc, type = "graph", layout = "spencer.circle", shownodesin = 3, vlabel = FALSE, node.pies=F)
# visualize node membership to different communities using node pies
plot(lc, type = "graph", shownodesin = 2, node.pies = F, vlabel = FALSE)
#visualize node community membership for the top-connected nodes using a community membership matrix
plot(lc, type = "members", nodes=unique(c(find.effects_TF$geneA, find.effects_TF$geneB)))
# display a summary of the results of the link communities algorithm
plot(lc, type = "summary", vlabel = FALSE)
# plot the dendrogram on its own with coloured community clusters
plot(lc, type = "dend")

# for nested communities
nested <- getAllNestedComm(lc)

getNestedHierarchies(lc, clusid = 1)
plot(lc, type = "graph", clusterids = c(9,11))



mc <- meta.communities(lc, hcmethod = "ward", deepSplit = 0)


# look for community centrality
cc <- getCommunityCentrality(lc)
head(sort(cc, decreasing = TRUE))

# community modularity and connectedness
# modularity of communities - the relative number of links within 
# the community versus links outside of the community
# and its inverse, community connectedness

cm <- getCommunityConnectedness(lc, conn = "modularity")
plot(lc, type = "commsumm", summary = "modularity")

# extract communities from dendogram at different height
lc2 <- newLinkCommsAt(lc, cutat = 0.4)

# extract nodes from communities
getNodesIn(lc, clusterids = c(4,5))

# find nodes shared by communities


