#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# The line below will be needed if you are running this script with python 2.
# Python 3 will ignore it.
from __future__ import print_function

# The following two lines will be needed in every python script:
from intermine.webservice import Service
# to get API key, go to https://yeastmine.yeastgenome.org/yeastmine/begin.do and create an account
# go to account details and create a new key

yeastmineAPItoken_file = open('/Users/Carolina/Documents/GitHub/DegreeProject/code/2020-02-05/yeastmineAPI.txt', 'r')
yeastmineAPItoken = yeastmineAPItoken_file.readline().rstrip()
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service", token = yeastmineAPItoken)

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "secondaryIdentifier", "symbol", "primaryIdentifier",
    "goAnnotation.ontologyTerm.identifier", "goAnnotation.ontologyTerm.name",
    "goAnnotation.ontologyTerm.parents.identifier",
    "goAnnotation.ontologyTerm.parents.name"
)

# You can edit the constraint values below
query.add_constraint("Gene", "IN", "genelist", code = "A")
query.add_constraint("organism.shortName", "=", "S. cerevisiae", code = "F")
query.add_constraint("status", "IS NULL", code = "C")
query.add_constraint("status", "=", "Active", code = "B")
query.add_constraint("goAnnotation.ontologyTerm.namespace", "=", "biological_process", code = "D")

# Your custom constraint logic is specified with the code below:
query.set_logic("(B or C) and F and A and D")
# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B")
print("gene", "symbol", "primaryIdentifier", \
    "ontologyTerm.identifier", "ontologyTerm.name", \
    "ontologyTerm.parents.identifier", \
    "ontologyTerm.parents.name", sep='\t')
for row in query.rows():
    print(row["secondaryIdentifier"], row["symbol"], row["primaryIdentifier"], \
        row["goAnnotation.ontologyTerm.identifier"], row["goAnnotation.ontologyTerm.name"], \
        row["goAnnotation.ontologyTerm.parents.identifier"], \
        row["goAnnotation.ontologyTerm.parents.name"], sep='\t')
