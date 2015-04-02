data(sotu)
head(sotu.tokens)
  
sotu.tokens = sotu.tokens[sotu.tokens$pos1 %in% c('N','M'),]

adjmat = wordWindowAdjacency(location=sotu.tokens$id,
                    term=sotu.tokens$lemma, 
                    context=sotu.tokens$aid,
                    window.size=10)

g = graph.adjacency(adjmat$adjmat, weighted=T)
V(g)$count = adjmat$termfreq

gs = getBackboneNetwork(g, alpha=0.05, max.vertices=100)
V(gs)$cluster = edge.betweenness.community(gs)$membership
gs = setNetworkAttributes(gs, size_attribute='count', cluster_attribute='cluster')
plot(gs)
