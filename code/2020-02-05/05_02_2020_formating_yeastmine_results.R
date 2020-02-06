library(data.table)

path <- "/Users/Carolina/Documents/GitHub/DegreeProject/" 
respath <- "/Users/Carolina/Documents/GitHub/DegreeProject/results/"
load(paste0(path,"results/2020-01-27/27_01_2020_linkedcommunities_TF.Rdata"))

genelist <- fread(paste0(path,"genelist.txt"), header = F)
cluster_table <- data.table(lc$nodeclusters)

im.yeast = initInterMine(listMines()["HumanMine"],token = Sys.getenv("yeastmine API"))


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
GOterms.count <- ympy.GO[, .(count = .N, description=unique(ontologyTerm.name)), by = ontologyTerm.identifier]
GOterms.count[order(-count)]


ympy.GO.cluster <- merge(ympy.GO, cluster_table, by.x="secondaryIdentifier", by.y="node", all=T, allow.cartesian=TRUE)

# find genes that have "transcription factor in the GO term"
ympy.GO.cluster[grepl("transcription factor |Transcription factor",goAnnotation.ontologyTerm.name, fixed = F)]




enrichResult = doEnrichment(
  im = im.yeast,
  ids = as.data.frame(genelist)$V1,
  widget = "go_enrichment_for_gene",
  filter = "biological_process",
  correction = "Benjamini Hochberg",
) 






# 
# 
# # results obtained at https://yeastmine.yeastgenome.org/yeastmine/bag.do?subtab=upload, 
# # using the name of the genes present in my network (gene list at /Users/Carolina/Documents/GitHub/DegreeProject/genelist.txt)
# # type = gene
# # organism = S. cerevisiae
# 
# 
# yeastmine_genefunctions <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastminesummary.csv")
# yeastmine_go <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastmine_go.csv")
# yeastmine_interactions.summary <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastmine_interactions_summary.tsv")
# yeastmine_interactions <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastmine_interactions.csv")
# 
# genefunctions.cluster <- merge(genefunctions, cluster_table, by.x="input", by.y="node", all.y=T)
# genefunctions.cluster[,organism.shortName:=NULL]
# genefunctions.cluster[,reason:=NULL]
# genefunctions.cluster[,matches:=NULL]
# 
# 
# genefunctions.cluster[cluster==2568]
