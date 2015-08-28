### Merge high similarity terms

#' Merge terms with a high mutual direction conditional probability
#' 
#' Merge terms that are likely to occur together. Specifically, terms of which the conditinal probability is higher than min.similarity in both directions. Or: P(A|B) >= min.similarity & P(B|A) >= min.similarity)
#' Not that this is not always a good thing. Unrelated terms that always occur together will also be merged. Whether this makes sense depends on the type of analysis.
#' 
#' @param m A sparse matrix in which columns are terms. Can be a DocumentTermMatrix class from the tm package
#' @param min.similarity The minimum conditional probability. The conditional probability of two terms in both directions needs to be higher than min.similarity for terms to be merged
#' @param max.label_length Terms that are merged together will be collapsed into a single label. To prevent very long labels, this is cut of from the [max.label_length] term.
#' @return a matrix (or document term matrix)
#' @export 
mergeTermClusters <- function(m, min.similarity=0.95, max.label_length=3){
  if('DocumentTermMatrix' %in% class(m)) {
    m = dtmToSparseMatrix(m) 
    return_dtm = T
  } else return_dtm = F
  message('Calculating term similarities')
  sim = termSimilarityNetwork(m, min.similarity)
  sim = decompose.graph(sim)
  N = length(sim)
  message('Merging ', N, ' word clusters')
  for(i in 1:N){
    if(i %% 1000 == 0) print(paste(i, '/', N))
    ids = as.numeric(V(sim[[i]])$name)
    label = collapseLabels(colnames(m)[ids], max.label_length)
    colnames(m)[ids] = label
  }
  
  m = collapseColumns(m)
  if(return_dtm) m = as.DocumentTermMatrix(m, weighting = weightTf)
  m
}

termSimilarityNetwork <- function(m, min.similarity){
  if('DocumentTermMatrix' %in% class(m)) m = dtmToSparseMatrix(m) 
  m@x[m@x > 0] = 1
  sim = crossprod(m) / colSums(m)
  
  sim@x = ifelse(sim@x >= min.similarity, 1, 0)
  sim = triu(triu(sim) + triu(t(sim))) # sum upper and lower triangle, and only take one
  sim = which(sim == 2, arr.ind = T) # select all word pairs where conditinal probability is higher than min.similarity in both ways (undirected because of lower tri)
  sim = sim[!sim[,1] == sim[,2],]
  graph.data.frame(sim, directed = F)
}

collapseLabels <- function(labels, max_length=3){
  if(length(labels) > max_length) {
    label = paste(paste(labels[1:max_length], collapse = '|'), paste('+', length(labels)-max_length, sep=''), sep='|')
  } else{
    label = paste(labels, collapse='|') 
  }
  label
}

collapseColumns <- function(m, as_mean=T){
  m = as(m, 'dgTMatrix')
  cnames = unique(colnames(m))
  j = match(colnames(m)[m@j+1], cnames)
  newm = spMatrix(nrow(m), length(cnames), m@i+1, j, m@x)
  newm = as(as(newm, 'dgCMatrix'), 'dgTMatrix')
  if(as_mean){
    cnames_count = table(colnames(m))
    cnames_count = as.numeric(cnames_count[match(cnames, names(cnames_count))])
    newm = t(t(newm) / cnames_count)
  }
  colnames(newm) = cnames
  rownames(newm) = rownames(m)
  newm
}

####

#' Reduce overlap of labels in igraph network
#' 
#' This function does not yet work properly. The idea is to use the wordlayout function of the wordcloud package to arrange nodes in order to prevent label overlap.
#' 
#' @param g A graph in the Igraph format, which has the g$layout attribute (set by using the layout. functions)
#' @param fontsize_multiplier The wordlayout function takes fontsize into account to determine overlap. As a temporary solution, increasing font size decreased label overlap (but messes up the network if done extremely)
#' @return a network layout, to be assigned to g$layout
#' @export 
reduceLabelOverlap <- function(g, fontsize_multiplier=5){
  layout_matrix = g$layout
  fontsize = V(g)$label.cex * fontsize_multiplier
  newlayout = wordlayout(layout_matrix[,1], layout_matrix[,2], V(g)$name, cex=fontsize)
  as.matrix(newlayout[,1:2])
}

#' Organize isolates next to network
#' 
#' Under development (if it proves usefull). Sometimes ignoring isolated words might seem bad, but having them placed randomly within and around the network is just plain annoying. This is an attempt to combine both sentiments.
#' 
#' @param g A graph in the Igraph format, which has the g$layout attribute (set by using the layout. functions)
#' @param fontsize_multiplier The wordlayout function takes fontsize into account to determine overlap. As a temporary solution, increasing font size decreased label overlap (but messes up the network if done extremely)
#' @return a network layout, to be assigned to g$layout
#' @export 
arrangeIsolates <- function(g){
  if(is.null(g$layout)) g$layout = layout.fruchterman.reingold(g)
  l = g$layout
  l[degree(g) == 0,1] = 0
  l[degree(g) == 0,2] = 0
  l = apply(l, 2, function(x) rescale(x, to=c(-100,100)))
  if(sum(degree(g)==0) == 0) return(l)
  
  isosizes = data.frame(id=1:sum(degree(g) == 0 ), size=V(g)$size[degree(g) == 0])
  isosizes = isosizes[order(-isosizes$size),]
  isosizes$y = 100 - cumsum(isosizes$size)
  isosizes = isosizes[order(isosizes$id),]
  
  isolayout = data.frame(x=rep(-130, nrow(isosizes)),y=isosizes$y)
  while(min(isolayout[,2]) < -100){
    isolayout[isolayout[,2] < -100,1] = isolayout[isolayout[,2] < -100,1] - 40
    isolayout[isolayout[,2] < -100,2] = isolayout[isolayout[,2] < -100,2] + 200
  }
  l[degree(g)==0,1] = isolayout[,1]
  l[degree(g)==0,2] = isolayout[,2]
  
  l = apply(l, 2, function(x) rescale(x, to=c(-100,100)))
  l
}

plotCluster <- function(g, cluster, redo_layout=F, ...){
  V(g)$cluster[is.na(V(g)$cluster)] = 0
  gs = delete.vertices(g, which(!V(g)$cluster == cluster))
  if(redo_layout) {gs$layout = layout.fruchterman.reingold(gs)
  }else gs$layout = gs$layout[V(g)$cluster == cluster,]
  plot(gs, ...)
  gs
}


