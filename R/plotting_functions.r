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

setVertexColors <- function(g, cluster){  
  V(g)$color = 'cadetblue1'
  if(!is.null(cluster)){
    pal = substr(rainbow(length(unique(cluster)), s=0.6,alpha=0.5), 1,7)
    duplicates = unique(cluster[duplicated(cluster)])
    cluster = match(cluster, duplicates) # re-index clusters, and setting isolates to NA
    V(g)$color[!is.na(cluster)] = pal[cluster[!is.na(cluster)]]
  }   
  V(g)$frame.color = V(g)$color
  g
}

setVertexAttributes <- function(g, size, cluster){
  g = setVertexColors(g, cluster)
  if(is.null(size)) {
    size = degree(g)
    message('No size attribute is given. Vertex size instead based on degree')
  }
  V(g)$size= rescale(size^0.4, to=c(2,15))
  V(g)$label.color = 'black'
  V(g)$label.cex = rescale(size, to=c(0.7,1.2))
  print('check')
  print(V(g)$name)
  V(g)$label = V(g)$name
  g
}

setEdgeAttributes <- function(g){
  E(g)$width = rescale(E(g)$weight, to=c(1,10))
  E(g)$arrow.size= E(g)$width/10
  E(g)$color='lightgrey'
  g
}


