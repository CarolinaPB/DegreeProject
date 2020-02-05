library(InterMineR)

# to get API key, go to https://yeastmine.yeastgenome.org/yeastmine/begin.do and create an account
# go to account details and create a new key
# create a file "~/.Renviron" and write yeastmineAPI = <token you just created>
# to access it Sys.getenv("yeastmine API")

# ran in the shell
# python /Users/Carolina/Documents/GitHub/DegreeProject/code/2020-02-05/yeastmineanalysis.py > ../../results/2020-02-05/yeastminepy/yeastmine_summar_py.txt








# results obtained at https://yeastmine.yeastgenome.org/yeastmine/bag.do?subtab=upload, 
# using the name of the genes present in my network (gene list at /Users/Carolina/Documents/GitHub/DegreeProject/genelist.txt)
# type = gene
# organism = S. cerevisiae


yeastmine_genefunctions <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastminesummary.csv")
yeastmine_go <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastmine_go.csv")
yeastmine_interactions.summary <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastmine_interactions_summary.tsv")
yeastmine_interactions <- fread("/Users/Carolina/Documents/GitHub/DegreeProject/results/2020-02-05/yeastmine_interactions.csv")

genefunctions.cluster <- merge(genefunctions, cluster_table, by.x="input", by.y="node", all.y=T)
genefunctions.cluster[,organism.shortName:=NULL]
genefunctions.cluster[,reason:=NULL]
genefunctions.cluster[,matches:=NULL]


genefunctions.cluster[cluster==2568]
