library(corpustools)
library(semnet)

data(sotu)
dtm = dtm.create(sotu.tokens$aid, sotu.tokens$lemma, sotu.tokens$freq, minfreq = 5)
dtm

g = chi2_wordassociations(dtm, fisher.p.thres=0.001)
test = fastgreedy.community(g)

V(g)$cluster = edge.betweenness.community(g, directed = F)$membership
g = setNetworkAttributes(g, size_attribute = V(g)$freq, cluster_attribute = V(g)$cluster)
plot(g, edge.arrow.size=0.0001)


data(sotu)
sotu.tokens = sotu.tokens[sotu.tokens$pos1 %in% c('M'),]
dtm = dtm.create(sotu.tokens$aid, sotu.tokens$lemma, sotu.tokens$freq, minfreq = 5)
dtm

d = chi2_wordassociations(dtm, chi.p.thres = 0.001, return.graph = F)
head(d)
min(d$odds_ratio)
head(d[order(-d$chi),])

d = chi2_wordassociations_alternative(dtm, p.thres = 0.001, return.graph = F)
min(d$odds_ratio)


g = chi2_wordassociations(dtm, fisher.p.thres=0.001)
V(g)$cluster = edge.betweenness.community(g, directed = F)$membership
g = setNetworkAttributes(g, size_attribute = V(g)$freq, cluster_attribute = V(g)$cluster)
plot(g, edge.arrow.size=0.0001)

library(semnet)
load('~/Dropbox/anp/compare_data.rdata')

d = chi2_wordassociations(dtm, p.thres = 0.001, ratio.thres = 1, return.graph = T)
deg = degree(g)
head(deg[order(-deg)], 50)




