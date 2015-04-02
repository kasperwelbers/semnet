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
