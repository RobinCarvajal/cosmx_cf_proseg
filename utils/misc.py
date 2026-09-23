import re

def grep(list, pattern):

    for item in list:
        if re.search(pattern, item, re.IGNORECASE):
            print(item)

## TEST
# from scripts.utils.misc import grep

# gene_names = list(comb.var_names)
# pattern = "^cd3"
# grep(gene_names, pattern)