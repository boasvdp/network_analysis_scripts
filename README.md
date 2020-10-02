# Network analysis scripts
Collection of network analysis scripts for ANI and/or SNP data

The main script at the moment is the `create_cytoscape_graph_from_distmat.R` script. This takes in a distance matrix (e.g. a SNP matrix or a matrix of %ANI values) and sets values below (ANI) or above (SNP) a certain threshold to 0. It then uses an igraph function to create a graph from this matrix and finally two functions defined at the bottom of the script, from https://github.com/idekerlab/cy-rest-R, to convert the igraph object to a `.json` that can be read using Cytoscape. 
