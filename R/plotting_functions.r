#' Set some default network attributes, including colors based on clustering, for pretty plotting
#' 
#' The purpose of this function is to create some default network attribute options to plot networks in a nice and insightfull way.
#' 
#' @param g A graph in the Igraph format.
#' @param size_attribute the name of the vertex attribute to be used to set the size of nodes
#' @param if TRUE, isolates are placed next to the network
#' @return a network in the Igraph format
#' @export 
setNetworkAttributes <- function(g, size_attribute=NULL, cluster_attribute=NULL){
  g = setVertexAttributes(g, size_attribute, cluster_attribute)
  g = setEdgeAttributes(g)  
  g$layout = layout.fruchterman.reingold(g)
  g
}

setVertexColors <- function(g, cluster_attribute){  
  if(!is.null(cluster_attribute)){
    cl = get.vertex.attribute(g, cluster_attribute)
    pal = colors(distinct=T)[grep('1', colors(distinct=T))] 
    pal = pal[grep('white|gray', pal, invert=T)] 
    pal = rep(pal, 5)
    
    duplicates = unique(cl[duplicated(cl)])
    cl = match(cl, duplicates) # re-index clusters, and setting isolates to NA
    V(g)$color = 'lightgrey'
    V(g)$color[!is.na(cl)] = pal[cl[!is.na(cl)]]
    V(g)$frame.color = V(g)$color
  } else {
    V(g)$color = 'cadetblue1'
    V(g)$frame.color = 'cadetblue1'
  }
  g
}

setVertexAttributes <- function(g, size_attribute, cluster_attribute){
  g = setVertexColors(g, cluster_attribute)
  proportion = if(!is.null(size_attribute)) get.vertex.attribute(g, size_attribute) else degree(g)
  V(g)$size= rescale(proportion^0.4, to=c(1,15))
  V(g)$label.color = 'black'
  V(g)$label.cex = rescale(proportion, to=c(0.7,1.2))
  V(g)$label = V(g)$name
  g
}

setEdgeAttributes <- function(g){
  E(g)$width = rescale(E(g)$weight, to=c(1,10))
  E(g)$arrow.size= E(g)$width/10
  E(g)$color='lightgrey'
  g
}

