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
yeastmineAPItoken_file = open('/Users/Carolina/Documents/GitHub/DegreeProject/code/2020-02-05/yeastmineAPI.txt', 'r')
yeastmineAPItoken = yeastmineAPItoken_file.readline().rstrip()
service = Service("https://yeastmine.yeastgenome.org/yeastmine/service", token = yeastmineAPItoken)

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view(
    "secondaryIdentifier", "chromosomeLocation.strand",
    "chromosomeLocation.start", "chromosomeLocation.end", "chromosome.primaryIdentifier"
)

# You can edit the constraint values below
query.add_constraint("Gene", "IN", "allgenes", code = "A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A")
print("gene", "strand", "start", "end", "chr", sep="\t")
for row in query.rows():
    print(row["secondaryIdentifier"], row["chromosomeLocation.strand"], \
        row["chromosomeLocation.start"], row["chromosomeLocation.end"], \
        row["chromosome.primaryIdentifier"], sep="\t")