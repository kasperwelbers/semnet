alphaForNVertices <- function(g, max.vertices){
  if(max.vertices > vcount(g)-2) return(1)
  a = get.data.frame(g)
  a = unique(data.frame(node=c(a$from,a$to), alpha=c(a$alpha,a$alpha)))
  a = a[order(a$alpha),]
  a = a[!duplicated(a$node),]
  max.alpha = a$alpha[max.vertices]
  if(max.alpha == a$alpha[max.vertices+1]) max.alpha = max.alpha - 0.000000001
  max.alpha
}

#' Extract the backbone of a network.
#' 
#' Based on the following paper: Serrano, M. Á., Boguñá, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the National Academy of Sciences, 106(16), 6483-6488.
#' 
#' @param g A graph in the `Igraph` format.
#' @param alpha The threshold for the alpha. Can be interpreted similar to a p value (see paper for clarrification). 
#' @param direction direction = 'none' can be used for both directed and undirected networks, and is (supposed to be) the disparity filter proposed in Serrano et al. (2009) is used. By setting to 'in' or 'out', the alpha is only calculated for out or in edges. This is an experimental use of the backbone extraction (so beware!) but it seems a logical application.  
#' @param delete.isolates If TRUE, vertices with degree 0 (i.e. no edges) are deleted.
#' @param max.vertices Optional. Set a maximum number of vertices for the network to be produced. The alpha is then automatically lowered to the point that only the given number of vertices remains connected (degree > 0). This can be usefull if the purpose is to make an interpretation friendly network. See e.g., http://jcom.sissa.it/archive/14/01/JCOM_1401_2015_A01 
#' @return A graph in the Igraph format
#' @export
getBackboneNetwork <- function(g, alpha=0.05, direction='none', delete.isolates=T, max.vertices=NULL){
  if(direction == 'none') E(g)$alpha = backbone.alpha(g)
  if(direction == 'in') E(g)$alpha = backbone.indegree.alpha(g)
  if(direction == 'out') E(g)$alpha = backbone.outdegree.alpha(g)
  g = delete.edges(g, which(E(g)$alpha >= alpha))
  if(!is.null(max.vertices) & ecount(g) > 0) {
    max.alpha = alphaForNVertices(g, max.vertices)
    if(max.alpha < alpha) message(paste('Using cutoff alpha', max.alpha, 'to keep N vertices under', max.vertices))
    g = delete.edges(g, which(E(g)$alpha >= max.alpha))
  }
  if(delete.isolates) g = delete.vertices(g, which(degree(g) == 0))
  if(ecount(g) == 0) {
    warning("No significant edges (backbone) remain!! Accept it (or lower the backbone_alpha)")
    return(g)
  }
  g
}

#' Calculate the alpha values that can be used to extract the backbone of a network.
#' 
#' Based on the following paper: Serrano, M. Á., Boguñá, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the National Academy of Sciences, 106(16), 6483-6488.
#' 
#' @param g A graph in the `Igraph` format.
#' @return A vector of alpha values, which matches the edges. Can thus easily be made an edge attribute: E(g)$alpha = backbone.alpha(g)
#' @export
backbone.alpha <- function(g){
  mat = get.adjacency(g, attr='weight')
  if(!is.directed(g)) mat[lower.tri(mat)] = 0 # prevents counting edges double in symmetric matrix (undirected graph)
  
  edgelist_ids = get.edgelist(g, names=F)
  alpha_ij = getAlpha(mat)[edgelist_ids] # alpha from the perspective of the 'from' node.
  alpha_ji = Matrix::t(getAlpha(Matrix::t(mat)))[edgelist_ids] # alpha from the perspective of the 'to' node.
  alpha_ij[alpha_ji < alpha_ij] = alpha_ji[alpha_ji < alpha_ij] # select lowest alpha, because an edge can be 'significant' from the perspective of both the 'from' and 'to' node. 
  alpha_ij
}

#' @export
getAlpha <- function(mat){
  weightsum = Matrix::rowSums(mat) + Matrix::colSums(mat)
  k = Matrix::rowSums(mat>0) + Matrix::colSums(mat>0)
  mat = mat / weightsum
  mat = (1 - mat)^(k-1)
  mat[is.na(mat)] = 1
  mat
}

#' Calculate the alpha values that can be used to extract the backbone of a network, for only the out.degree
#' 
#' Based on the following paper: Serrano, M. Á., Boguñá, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the National Academy of Sciences, 106(16), 6483-6488.
#' 
#' @param g A graph in the `Igraph` format.
#' @return A vector of alpha values, which matches the edges. Can thus easily be made an edge attribute: E(g)$alpha = backbone.alpha(g)
#' @export
backbone.outdegree.alpha <- function(g){
  mat = get.adjacency(g, attr='weight')
  weightsum = Matrix::rowSums(mat)
  k = Matrix::rowSums(mat > 0)
  edgelist_ids = get.edgelist(g, names=F)
  getAlpha(mat, weightsum, k)[edgelist_ids]
}

#' Calculate the alpha values that can be used to extract the backbone of a network, for only the in.degree
#' 
#' Based on the following paper: Serrano, M. Á., Boguñá, M., & Vespignani, A. (2009). Extracting the multiscale backbone of complex weighted networks. Proceedings of the National Academy of Sciences, 106(16), 6483-6488.
#' 
#' @param g A graph in the `Igraph` format.
#' @return A vector of alpha values, which matches the edges. Can thus easily be made an edge attribute: E(g)$alpha = backbone.alpha(g)
#' @export
backbone.indegree.alpha <- function(g){
  mat = get.adjacency(g, attr='weight')
  weightsum = Matrix::colSums(mat)
  k = Matrix::colSums(mat > 0)
  edgelist_ids = get.edgelist(g, names=F)
  t(getAlpha(t(mat), weightsum, k))[edgelist_ids]
}
