data(sotu)
head(sotu.tokens)
  
sotu.tokens = sotu.tokens[sotu.tokens$pos1 %in% c('N','M'),]

adj = wordWindowAdjacency(location=sotu.tokens$id,
                    term=sotu.tokens$lemma, 
                    context=sotu.tokens$aid,
                    window.size=10)

g = graph.adjacency(adj$adj, weighted=T)
V(g)$count = adj$termfreq

gs = getBackboneNetwork(g, alpha=0.05, max.vertices=100)
V(gs)$cluster = edge.betweenness.community(gs)$membership
gs = setNetworkAttributes(gs, size_attribute='count', cluster_attribute='cluster')
plot(gs)

## get co-occurence of words per context

# first select specific words (using all words will result in a huuuuuuuuge data.frame)
sotu.tokens2 = sotu.tokens[grepl('health|care|reform|plan', sotu.tokens$lemma),]
adj = wordWindowAdjacency(location=sotu.tokens2$id,
                             term=sotu.tokens2$lemma, 
                             context=sotu.tokens2$aid,
                             window.size=10,
                             output.per.context=T)
head(adj$adj, 20)

