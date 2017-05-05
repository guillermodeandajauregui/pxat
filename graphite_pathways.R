#use graphite + igraph to analyze pathways 

library(graphite)
library(igraph)
library(org.Hs.eg.db)

reactome_graphites = graphite::pathways(species = "hsapiens", database = "reactome") #makes list of pathways
reactome_graphites = lapply(reactome_graphites, convertIdentifiers, "symbol")
reactome_graphites = reactome_graphites[1:10]
reactome_graphites_graphs = lapply(X = reactome_graphites, FUN = graphite::pathwayGraph)
reactome_graphites_eyegraphs = lapply(reactome_graphites_graphs, FUN = igraph::graph_from_graphnel)

eyegraph_union = igraph::union(reactome_graphites_eyegraphs[[1]], reactome_graphites_eyegraphs[[2]])
for(i in 3:length(reactome_graphites_eyegraphs)){
  eyegraph_union = igraph::union(eyegraph_union, reactome_graphites_eyegraphs[[i]])
}
plot(eyegraph_union)
