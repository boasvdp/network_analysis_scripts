## This code reads a distance matrix (e.g. created by snp-dists)
## This distmat is then read as dataframe with row and column names
## Values in distmat matrix are numbers of differences between strains

# load libraries
library(jsonlite)
library(igraph)

# Read distmat file 
distmat <- read.delim("~/projects/thomas_sysrev/fastANI_1719_distmat.tsv", row.names=1)

# filter distmat based on your criteria
# EDGES WITH VALUE ZERO WILL BE REMOVED
# IF YOU USE SNPS AND SOME STRAINS ARE 0 SNPS APART: SET THESE TO 1 TO PREVENT FILTERING
# Set 0's to 1's to prevent removal
#distmat[distmat<1]=1
# Set SNP distances of more than 10 to 0, to get them removed
distmat[distmat<95]=0

# Convert dataframe to matrix
dm <- data.matrix(distmat, rownames.force = NA)

# Make a graph from the matrix. Function is custom, create this with the code at the bottom!
graph <- graph_from_adjacency_matrix(dm, mode = "undirected", weight = TRUE, diag = FALSE, add.colnames = NA, add.rownames = NULL)

# Convert to cytoscape object. Function is custom, create this with the code at the bottom!
g <- toCytoscape(graph)

# write g to a local .json file
write(g,"~/projects/thomas_sysrev/graph_1719_suis_ANI96.json")

##################################################################################
## toCytoscape.R
toCytoscape <- function (igraphobj) {
  # Extract graph attributes
  graph_attr = graph.attributes(igraphobj)
  
  # Extract nodes
  node_count = length(V(igraphobj))
  if('name' %in% list.vertex.attributes(igraphobj)) {
    V(igraphobj)$id <- V(igraphobj)$name
  } else {
    V(igraphobj)$id <- as.character(c(1:node_count))
  }
  
  nodes <- V(igraphobj)
  v_attr = vertex.attributes(igraphobj)
  v_names = list.vertex.attributes(igraphobj)
  
  nds <- array(0, dim=c(node_count))
  for(i in 1:node_count) {
    if(i %% 1000 == 0) {
      print(i)
    }
    nds[[i]] = list(data = mapAttributes(v_names, v_attr, i))
  }
  
  edges <- get.edgelist(igraphobj)
  edge_count = ecount(igraphobj)
  e_attr <- edge.attributes(igraphobj)
  e_names = list.edge.attributes(igraphobj)
  
  attr_exists = FALSE
  e_names_len = 0
  if(identical(e_names, character(0)) == FALSE) {
    attr_exists = TRUE
    e_names_len = length(e_names)
  }
  e_names_len <- length(e_names)
  
  eds <- array(0, dim=c(edge_count))
  for(i in 1:edge_count) {
    st = list(source=toString(edges[i,1]), target=toString(edges[i,2]))
    
    # Extract attributes
    if(attr_exists) {
      eds[[i]] = list(data=c(st, mapAttributes(e_names, e_attr, i)))
    } else {
      eds[[i]] = list(data=st)
    }
    
    if(i %% 1000 == 0) {
      print(i)
    }
  }
  
  el = list(nodes=nds, edges=eds)
  
  x <- list(data = graph_attr, elements = el)
  #print("Done.")
  return (toJSON(x, auto_unbox=TRUE, pretty=TRUE))
}
###############################################################################

###############################################################################
## mapAttributes function
mapAttributes <- function(attr.names, all.attr, i) {
  attr = list()
  cur.attr.names = attr.names
  attr.names.length = length(attr.names)
  
  for(j in 1:attr.names.length) {
    if(is.na(all.attr[[j]][i]) == FALSE) {
      #       attr[j] = all.attr[[j]][i]
      attr <- c(attr, all.attr[[j]][i])
    } else {
      cur.attr.names <- cur.attr.names[cur.attr.names != attr.names[j]]
    }
  }
  names(attr) = cur.attr.names
  return (attr)
}
################################################################################
