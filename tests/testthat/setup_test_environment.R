## For some reason, some of the packages in Suggests are not found while 
## running the unit tests during R CMD check. This seems to help.
requireNamespace("wordspace")
requireNamespace("tidyr")
requireNamespace("Gviz")
