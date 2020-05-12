# run with python3
import os.path
from os import path
os.chdir("/Users/Carolina/Documents/GitHub/DegreeProject/summarizing/")
if not path.isfile("results/essentialgenes.txt"):
  # The following lines will be needed in every python script:
  from intermine.webservice import Service
  yeastmineAPItoken_file = open('data/yeastmineAPI.txt', 'r')
  yeastmineAPItoken = yeastmineAPItoken_file.readline().rstrip()
  service = Service("https://yeastmine.yeastgenome.org/yeastmine/service", token = yeastmineAPItoken)

  # Get a new query on the class (table) you will be querying:
  query = service.new_query("Gene")

  # The view specifies the output columns
  query.add_view("secondaryIdentifier", "symbol", "name", "phenotypeSummary")

  # You can edit the constraint values below
  query.add_constraint("Gene", "IN", "allgenes", code = "A")

  terms = "gene", "symbol", "gene.name", "phenotype.summary"

  terms_query = ["secondaryIdentifier", "symbol", "name", "phenotypeSummary"]
  print("saving file")
  with open("results/essentialgenes.txt", "w") as file:
    # write headers
    for term in terms[:-1]:
      file.write(term)
      file.write("\t")
    else:
      file.write(terms[-1])
    file.write("\n")
    #write content
    for row in query.rows():
      for t in terms_query[:-1]:
        if row[t] != None:
          file.write(str(row[t]))
          file.write("\t")
        else:
          file.write("NA")
          file.write("\t")
      if row[terms_query[-1]] != None:
        file.write(row[terms_query[-1]])
      else:
        file.write("NA")
      file.write("\n")
else:
  print("file exists")
