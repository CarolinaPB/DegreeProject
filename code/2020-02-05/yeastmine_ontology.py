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

yeastmineAPItoken_file = open('yeastmineAPI.txt', 'r')
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
query.add_constraint("goAnnotation.ontologyTerm.parents.identifier", "ONE OF", ["GO:1901362", "GO:0018130", "GO:0019438", "GO:0034654", "GO:0044281", "GO:0006629", "GO:1903506", "GO:2001141", "GO:0019219", "GO:0051252", "GO:1901575", "GO:0006355"], code = "B")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B")
print("secondaryIdentifier", "symbol", "primaryIdentifier", \
    "goAnnotation.ontologyTerm.identifier", "goAnnotation.ontologyTerm.name", \
    "goAnnotation.ontologyTerm.parents.identifier", \
    "goAnnotation.ontologyTerm.parents.name", sep='\t')
for row in query.rows():
    print(row["secondaryIdentifier"], row["symbol"], row["primaryIdentifier"], \
        row["goAnnotation.ontologyTerm.identifier"], row["goAnnotation.ontologyTerm.name"], \
        row["goAnnotation.ontologyTerm.parents.identifier"], \
        row["goAnnotation.ontologyTerm.parents.name"], sep='\t')
