#' Visualize a semnet network
#' 
#' plot_semnet is a wrapper for the plot.igraph() function optimized for plotting a network with vertex labels. 
#' 
#' Before plotting the network, the setNetworkAttributes() function is used to set pretty defaults for plotting. Optionally, reduce_labeloverlap can be used to prevent labeloverlap (as much as possible).
#'
#' @param g A network in the igraph format. Specifically designed for the output of coOccurenceNetwork() and windowedCoOccurenceNetwork()
#' @param vertexsize_attr a character string indicating a vertex attribute that represents size. Default is 'freq', which is created in the coOccurenceNetwork functions to indicate the number of times a word occured. 
#' @param vertexcolor_attr a character string indicating a vertex attribute that represents color. The attribute can also be a numeric value (e.g., a cluster membership) in which case colors are assigned to numbers. 
#' @param labelsize_coef a coefficient for increasing or decreasing the size of the vertexlabel.
#' @param labeldist_coef a coefficient that roughly determines the minimal distance between vertex labels, based on the size of labels. Only used if reduce_labeloverlap is TRUE.
#' @param reduce_labeloverlap if TRUE, an algorithm is used to reduce overlap as best as possible.
#' @param redo_layout If TRUE, a new layout will be calculated using layout_with_fr(). If g does not have a layout attribute (g$layout), a new layout is automatically calculated.
#' @param return_graph if TRUE, plot_semnet() also returns the graph object with the attributes and layout as shown in the plot. 
#' @param ... additional arguments are passed on to plot.igraph()
#'
#' @return Plots a network, and returns the network object if return_graph is TRUE.
#' @export 
plot_semnet <- function(g, weight_attr='weight', minweight=NA, vertexsize_attr='freq', vertexcolor_attr=NA, labelsize_coef=1, labeldist_coef=1.5, reduce_labeloverlap=T, redo_layout=F, return_graph=T, layout=layout_with_fr, ...){
  if(redo_layout | is.null(g$layout)) g$layout = layout(g)
  
  E(g)$weight = get.edge.attribute(g, weight_attr)
  if(!is.na(minweight)) g = delete.edges(g, which(E(g)$weight < minweight))
  
  g = setNetworkAttributes(g, vertexsize_attr, vertexcolor_attr)
  V(g)$label.cex = labelsize_coef * V(g)$label.cex
  
  if(reduce_labeloverlap){
    g = reduceLabelOverlap(g, labeldist_coef, cex_from_device = T) 
  }
  g = plotArgsToAttributes(g, args=list(...))

  plot.igraph(g, ...)
  
  if(return_graph) g
}


plotArgsToAttributes <- function(g, args){
  if(length(args) == 0) return(g)
  for(i in 1:length(args)){
    name = names(args)[i]
    if(!grepl('vertex\\.|edge\\.', name)) g = set.graph.attribute(g, name, value=args[[i]])
    if(grepl('vertex\\.', name)) g = set.vertex.attribute(g, gsub('vertex\\.', '', name), value=args[[i]])
    if(grepl('edge\\.', name)) g = set.edge.attribute(g, gsub('edge\\.', '', name), value=args[[i]])
  }
  g
}

#' Set some default network attributes for pretty plotting
#' 
#' The purpose of this function is to create some default network attribute options to plot networks in a nice and insightfull way.
#' 
#' @param g A graph in the Igraph format.
#' @param size_attribute the name of the vertex attribute to be used to set the size of nodes
#' @param if TRUE, isolates are placed next to the network
#' @return a network in the Igraph format
#' @export 
setNetworkAttributes <- function(g, size_attribute='freq', color_attribute=NA, redo_layout=F){
  g = setVertexAttributes(g, size_attribute, color_attribute)
  g = setEdgeAttributes(g)  
  if(is.null(g$layout)) g$layout = layout_with_fr(g)
  g
}

setVertexColors <- function(g, color){  
  if(!is.null(color)){
    if(class(color) == 'numeric'){
      pal = substr(rainbow(length(unique(color)), s=0.6,alpha=0.5), 1,7)
      duplicates = unique(color[duplicated(color)])
      color = match(color, duplicates) # re-index colors, and setting isolates to NA
      V(g)$color[!is.na(color)] = pal[color[!is.na(color)]]
    } else {
      V(g)$color = color
    }
  } else {
    V(g)$color = 'cadetblue1' 
  }
  V(g)$frame.color = V(g)$color
  g
}

setVertexAttributes <- function(g, size, color){
  vattrs = names(vertex.attributes(g))
  if(is.na(color) | !color %in% vattrs) {
    color = fastgreedy.community(as.undirected(g))$membership
    message('No (valid) color attribute given. Vertex color now based on undirected fastgreedy.community() clustering')
  } else {
    color = get.vertex.attribute(g, color)
  }
  g = setVertexColors(g, color)
  
  if(is.na(size) | !size %in% vattrs) {
    size = rep(1, vcount(g))
    message('No (valid) size attribute given. Vertex size now set to 1')
  } else {
    size = get.vertex.attribute(g, size)
  }
  
  V(g)$size= rescale_var(size^0.4, 2, 15)
  V(g)$label.color = 'black'
  V(g)$label.cex = rescale_var(size, 0.85, 1.25)
  V(g)$label = V(g)$name
  g
}


setEdgeAttributes <- function(g){
  E(g)$width = rescale_var(E(g)$weight, 1, 10)
  E(g)$arrow.size= E(g)$width/10
  E(g)$color='lightgrey'
  g
}


rescale_var <- function(x, new_min=0, new_max=1, x_min=min(x), x_max=max(x)){
  if(x_min == x_max) return(x)
  x = (x - x_min) / (x_max - x_min) # normalize
  x = x * (new_max-new_min)
  return(x + new_min)
}

reduceLabelOverlap <- function(g, labeldist_coef=1.5, cex_from_device=F, label.attr='label', labelsize.attr='label.cex'){
  layout_matrix = layout.norm(g$layout)
  
  vnames = names(vertex.attributes(g))
  if(!label.attr %in% vnames) {
    stop('"', label.attr, '" is not a valid vertex attribute)')
  } else {
    label = get.vertex.attribute(g, label.attr)
  }
  
  if(labelsize.attr %in% vnames){
    label.cex = get.vertex.attribute(g, labelsize.attr)
  } else {
    message('"', labelsize.attr, '" is not a valid vertex attribute). Labelsize is set to 1')
    label.cex = 1
  }
  
  ord = order(-degree(g)) ## reorder so that words with the least connections are rearranged first
  layout_matrix = layout_matrix[ord,]
  label.cex = label.cex[ord]
  label = label[ord]
 
  plot(range(layout_matrix), axes = F, frame.plot = F, xlab='', ylab='', type='n')
  newlayout = wordlayout(layout_matrix[,1], layout_matrix[,2], label, cex=label.cex*(labeldist_coef), rstep=0.05)
  
  ## calculate new cex based on percentual difference old and new word width
  #oldwidth = strwidth(label, cex=label.cex*(labeldist_coef))
  #shrinkcoef = newlayout[,'width'] / oldwidth
  #newlayout = cbind(newlayout, newcex=label.cex*shrinkcoef)
  
  newlayout = newlayout[match(1:nrow(newlayout), ord),] # return original order
    
  #g = set.vertex.attribute(g, labelsize.attr, value=newlayout[,'newcex'])
  g$layout = as.matrix(newlayout[,1:2])
  g
}

