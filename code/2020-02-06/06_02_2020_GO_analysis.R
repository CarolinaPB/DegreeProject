library(data.table)
library(linkcomm) # for the plot

path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
load(paste0(path,"results/2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))

genelist <- fread(paste0(path,"genelist.txt"))
cluster_table <- data.table(lc$nodeclusters)

im.yeast = initInterMine(listMines()["HumanMine"],token = Sys.getenv("yeastmine API"))

find.effects_TF <- fread(paste0(respath, "2020-01-27/findeffects_TF.gz"))

# to get API key, go to https://yeastmine.yeastgenome.org/yeastmine/begin.do and create an account
# go to account details and create a new key
# create a file "~/.Renviron" and write yeastmineAPI = <token you just created>
# to access it Sys.getenv("yeastmine API")

#### create the yeastmine files in Python ####

# to create the summary results table
# python /Users/Carolina/Documents/GitHub/DegreeProject/code/2020-02-05/yeastmineanalysis.py > /Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastminepy/yeastmine_summary_py.txt
ympy.summary <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastminepy/yeastmine_summary_py.txt")

# to create the GO table
# python /Users/Carolina/Documents/GitHub/DegreeProject/code/2020-02-05/yeastmine_ontology.py > /Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastminepy/yeastmine_ontology_py.txt
ympy.GO <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastminepy/yeastmine_ontology_py.txt")

# to create enrichment table
# python /Users/Carolina/Documents/GitHub/DegreeProject/code/2020-02-06/yeastminepy_goenrichment.py > /Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-06/enrichment.txt
ympy.enrichment <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-06/enrichment.txt")
############

## count number of times each ontology term appears
GOterms.count <- ympy.GO[, .(description=unique(ontologyTerm.name), count = .N), by = ontologyTerm.identifier]
GOterms.count[order(-count)]

# add cluster number to the table with gene and GO term
ympy.GO.cluster <- merge(ympy.GO, cluster_table, by.x="gene", by.y="node", all=T, allow.cartesian=TRUE)

# find genes that have "transcription factor" in the GO term
GO_trans_factor <- ympy.GO.cluster[grepl("transcription factor",ontologyTerm.name, fixed = F)]
unique(GO_trans_factor[,.(gene, ontologyTerm.identifier, ontologyTerm.name)])

# number of times each geneA points to another gene
numlinks_from_gene <- find.effects_TF[, .(count=.N), by=geneA]

# add each time that gene points to another
links_perGO_pergene <- merge(unique(ympy.GO[,.(gene, ontologyTerm.name, ontologyTerm.identifier)]), numlinks_from_gene, by.x="gene", by.y="geneA", all=T)

# table with GO identifier, description and number of times it belongs to genes pointing to other genes
links_perGO <- links_perGO_pergene[,.(gocount=sum(count, na.rm = T)), by=c("ontologyTerm.name", "ontologyTerm.identifier")]

# find the GOs that have transcription in the name
links_perGO_transfactor <- links_perGO[grep(pattern = "transcription factor", ontologyTerm.name)]
links_perGO_transfactor <- links_perGO_transfactor[order(-gocount)]

# barplot of number of links coming out of genes that have GO term with "transcription" factor in it
#pdf("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-06/link_per_transcfactorGOterm.pdf")
bp <- barplot(links_perGO_transfactor$gocount, names.arg = links_perGO_transfactor$ontologyTerm.identifier, 
        ylab = "# links", xlab ="GO" , ylim = c(0,200), main="# links going out of genes with GO term \n only GO terms including 'transcription factor'")
    
text(bp,links_perGO_transfactor$gocount,links_perGO_transfactor$gocount,cex=0.8, pos=3)    
#dev.off()
links_perGO_transfactor


# table with geneA, geneB and the GO term for gene A
find.effects_TF_GO_geneA <- merge(unique(find.effects_TF[,.(geneA, geneB)]), unique(ympy.GO[,.(gene, ontologyTerm.identifier)]), all=T, by.x="geneA", by.y="gene", allow.cartesian=TRUE)
# table with geneA, geneB, GO term and cluster number for the gene pair
find.effects_TF_GO_geneA.cluster <- merge(find.effects_TF_GO_geneA, data.table(lc$edges), by.x=c("geneA", "geneB"), by.y=c("node1", "node2"), allow.cartesian = T, all=T)

find.effects_TF_GO_geneA.cluster[!is.na(cluster)] 


ids_to_plot <- unique(find.effects_TF_GO_geneA.cluster[ontologyTerm.identifier %in% links_perGO_transfactor$ontologyTerm.identifier]$cluster)

plot(lc, type = "graph", node.pies = F, vlabel=F, vshape="circle" , vsize=3, ewidth=1, clusterids= ids_to_plot,arrow.size=0.01)


load(paste0(respath, "2020-01-29/communitymodularity.Rdata"))

con_mod <- data.table(cluster=1:length(cm),connectedness=cconnect, modularity=cm)

con_mod_numnodes <- merge(con_mod,data.table(cluster=as.integer(names(lc$clustsizes)), numnodes=lc$clustsizes))
con_mod_numnodes <- con_mod_numnodes[order(-modularity)]
topmod_clusters <- con_mod_numnodes$cluster[1:10]
topmod_clusters.name <- paste("cluster_",topmod_clusters, sep="")

length(find.effects_TF_GO_geneA.cluster[cluster %in% topmod_clusters]$ontologyTerm.identifier)

find.effects_TF_GO_geneA.cluster[cluster==2358]

ympy.GO.cluster.sub <- unique(ympy.GO.cluster[,.(gene, ontologyTerm.identifier, ontologyTerm.name, cluster)])



#### Analyse cluster with higher modularity
ympy.GO.cluster.sub[cluster==2358]


plot(lc, type = "graph", node.pies = F, vlabel=T, vshape="circle" , vsize=3, ewidth=1, clusterids= 2358,arrow.size=0.01)


# YEL055C affects YEL056W and YEL067C
# YEL056W seems to be involved in chromatine stuff 
# YEL055C is involved in ribosome genesis and synthesis of RNA
# YEL067C function unknown

# YEL036C membrane, golgi  and stuff
# 

# makes more sense that YEL055C points at YEL056W and YEL067C than YEL036C pointing at those

#### Analyse cluster with second highes modularity
ympy.GO.cluster.sub[cluster==2568]


plot(lc, type = "graph", node.pies = F, vlabel=T, vshape="circle" , vsize=3, ewidth=1, clusterids= 2568,arrow.size=0.01)

