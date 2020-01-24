demo(topic = "linkcomm", package = "linkcomm")

lc <- getLinkCommunities(find.effects_TF[,.(geneA, geneB)], directed = T)
print(lc)

plot(lc, type = "graph", layout = layout.fruchterman.reingold, vlabel = FALSE)

# The following displays only the nodes that belong to 3 or more communities:
plot(lc, type = "graph", layout = "spencer.circle", shownodesin = 3)
plot(lc, type = "graph", shownodesin = 2, node.pies = TRUE)
plot(lc, type = "members")
plot(lc, type = "summary")
plot(lc, type = "dend")

# for nested communities
getAllNestedComm(lc)

getNestedHierarchies(lc, clusid = 9)
plot(lc, type = "graph", clusterids = c(9,11))

cr <- getClusterRelatedness(lc, hcmethod = "ward")
cutDendrogramAt(cr, cutat = 1.2)

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