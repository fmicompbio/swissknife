print(requireNamespace("wordspace"))
print(requireNamespace("tidyr"))
print(requireNamespace("Gviz"))
print(requireNamespace("SingleR"))

install.packages("sessioninfo")
pkgs <- installed.packages()[, "Package"]
print(sessioninfo::session_info(pkgs, include_base = TRUE))