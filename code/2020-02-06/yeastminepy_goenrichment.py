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

lm=service.list_manager()
l=lm.get_list(name="genelist")

r=l.calculate_enrichment(widget="go_enrichment_for_gene",maxp=0.05)

print("GO", "description", "pval", "nmatches", sep="\t")
for i in r:
    print(i.identifier,i.description,i.p_value, i.matches, sep="\t")
