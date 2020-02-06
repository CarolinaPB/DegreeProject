# Intermine snippets


# To view the names of all the available templates you can view the templates dictionary and iterate through it.
service = Service("www.flymine.org/flymine/service")
templates=service.templates
for name in templates.keys():
    print(name)
