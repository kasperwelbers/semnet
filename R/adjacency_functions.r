#' Create an adjacency graph from a document term matrix
#' 
#' @param dtm a document term matrix, either or not in the dtm format from the `tm` package
#' @param measure the measure to calcualte adjacency. Currently supports cosine and conditional probability
#' @return A graph in the Igraph format in which edges represent the adjacency of terms
#' @export
adjacencyGraph <- function(dtm, measure='cosine'){
  if('DocumentTermMatrix' %in% class(dtm)) dtm = dtmToSparseMatrix(dtm)
  if(measure == 'cosine') g = graph.adjacency(getCosine(dtm), mode='undirected', weighted=T)
  if(measure == 'conprob') g = graph.adjacency(getConditionalProbability(dtm), mode='directed', weighted=T)
  V(g)$occurence = col_sums(dtm) / nrow(dtm)
  g
}

dtmToSparseMatrix <- function(dtm){
  sm = spMatrix(nrow(dtm), ncol(dtm), dtm$i, dtm$j, dtm$v)
  rownames(sm) = rownames(dtm)
  colnames(sm) = colnames(dtm)
  sm
}

getConditionalProbability <- function(mat){
  mat@x[mat@x > 0] = 1
  mat = Matrix::crossprod(mat)/colSums(mat) # conditional probability of words, based on the co-occurence of words in the same document AND same topic
  Matrix::diag(mat) = 0 # ignore loops
  mat[is.na(mat)] = 0
  as.matrix(mat) ## Without changing it to regular matrix, I get an error for large matrices when converting it to an Igraph graph. (Though this is not the case if run line-by-line instead of using this function, which suggests some sort of local/global issue)
}

getCosine <- function(mat){
  mat = Matrix::crossprod(mat)
  mat = mat/Matrix::crossprod(t(Matrix::diag(sqrt(mat))))
  Matrix::diag(mat) = 0 # ignore loops
  mat[is.na(mat)] = 0
  as.matrix(mat) ## Without changing it to regular matrix, I get an error for large matrices when converting it to an Igraph graph. (Though this is not the case if run line-by-line instead of using this function, which suggests some sort of local/global issue)
}



##### windowed adjacency functions #####

stretchLocation <- function(location, context, window.size){
  ## (location and context need to be sorted on order(location, context))
  newcontext = which(!duplicated(context))
  context.max = location[newcontext-1] + (window.size*2)
  multiplier_scores = cumsum(c(0,context.max))
  multiplier_vector = rep(NA, length(context))
  multiplier_vector[newcontext] = multiplier_scores
  multiplier_vector = na.locf(multiplier_vector)
  return(location + multiplier_vector)
}

locationMatrix <- function(i, j, shifts=0, count.double=F){
  mat = spMatrix(max(i), max(j))
  for(shift in shifts){
    i_shift = i + shift
    select = i_shift > 0 & i_shift <= max(i)
    mat = mat + spMatrix(nrow=max(i), ncol=max(j), i=i_shift[select], j=j[select], rep(1, sum(select))) 
  }
  mat = mat[i,]
  if(count.double==F) mat[mat>0] = 1
  mat
}

locationWindowMatrix <- function(location, term, context, window.size=3, two.sided=T, count.double=F){
  ord = order(context, location)
  location = location[ord,]
  term = term[ord,]
  context = context[ord,]
  
  location = stretchLocation(location,context,window.size=3)
  shifts = if(two.sided) -window.size:window.size else 0:window.size
  terms = unique(term)
  term_index = match(term, terms)
  
  location.mat = locationMatrix(location, term_index, 0)
  window.mat = locationMatrix(location, term_index, shifts, count.double)
  colnames(location.mat) = colnames(window.mat) = terms
  rownames(location.mat) = rownames(window.mat) = context
  list(location.mat=location.mat, window.mat=window.mat)
}








