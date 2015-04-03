data(sotu)
head(sotu.tokens)
  
sotu.tokens = sotu.tokens[sotu.tokens$pos1 %in% c('N','M'),]

adj = wordWindowAdjacency(location=sotu.tokens$id,
                    term=sotu.tokens$lemma, 
                    context=sotu.tokens$aid,
                    window.size=10)

g = graph.adjacency(adj, weighted=T)

gs = getBackboneNetwork(g, alpha=0.05, max.vertices=100)
V(gs)$cluster = edge.betweenness.community(gs)$membership
gs = setNetworkAttributes(gs, cluster_attribute='cluster')
plot(gs)


## get co-occurence of words per context

# first select specific words (using all words will result in a huuuuuuuuge data.frame)
sotu.tokens2 = sotu.tokens[grepl('health|care|reform|plan', sotu.tokens$lemma),]
adj = wordWindowAdjacency(location=sotu.tokens2$id,
                             term=sotu.tokens2$lemma, 
                             context=sotu.tokens2$aid,
                             window.size=10,
                             output.per.context=T)
head(adj, 20)


## merge adjacency matrices. (solution for big data limitations)

batch1 = sotu.tokens[1:1000,]
batch2 = sotu.tokens[1001:2000,]
batch3 = sotu.tokens[2001:3000,]

adj1 = wordWindowAdjacency(batch1$id, batch1$lemma, batch1$aid, window.size=10)
adj2 = wordWindowAdjacency(batch2$id, batch2$lemma, batch2$aid, window.size=10)
adj3 = wordWindowAdjacency(batch3$id, batch3$lemma, batch3$aid, window.size=10)

allmatrices = list(adj1,adj2,adj3)
adj = mergeMatrices(allmatrices)

