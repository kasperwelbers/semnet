

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
  ## (location and context need to be sorted on order(context,location))
  newcontext = which(!duplicated(context))
  context.max = location[newcontext-1] + (window.size*2)
  multiplier_scores = cumsum(c(0,context.max))
  multiplier_vector = rep(NA, length(context))
  multiplier_vector[newcontext] = multiplier_scores
  multiplier_vector = na.locf(multiplier_vector)
  return(location + multiplier_vector)
}

locationMatrix <- function(i, j, shifts=0, count.once=T){
  mat = spMatrix(max(i), max(j))
  for(shift in shifts){
    i_shift = i + shift
    select = i_shift > 0 & i_shift <= max(i)
    mat = mat + spMatrix(nrow=max(i), ncol=max(j), i=i_shift[select], j=j[select], rep(1, sum(select))) 
  }
  mat = mat[i,]
  if(count.once) mat[mat>0] = 1
  mat
}

#' A sliding window approach to calculate the co-occurence of words
#' 
#' @param location An integer vector giving the position of terms in a given context (e.g., document, paragraph, sentence) 
#' @param term A character vector giving the terms
#' @param context A vector giving the context in which terms occur (e.g., document, paragraph, sentence)
#' @param window.size The distance within which words should occur from each other to be counted as a co-occurence.
#' @param two.sided Logical. If false, it is only counted how often a word occured `after` another word within the given window size
#' @param count.once Logical. A word can occur multiple times within the same window. If count.once is true, this is only counted as a single occurence (suggested).
#' @return A graph in the Igraph format in which edges represent the adjacency of terms
#' @export
wordWindowAdjacency <- function(location, term, context, window.size=3, two.sided=T, count.once=T){
  ord = order(context, location)
  location = location[ord]
  term = term[ord]
  context = context[ord]
  
  location = stretchLocation(location,context,window.size=3)
  shifts = if(two.sided) -window.size:window.size else 0:window.size
  terms = unique(term)
  term_index = match(term, terms)
  
  location.mat = locationMatrix(location, term_index, 0)
  window.mat = locationMatrix(location, term_index, shifts[!shifts == 0], count.once)
  colnames(location.mat) = colnames(window.mat) = terms
  
  calculateAdjacency(location.mat, window.mat)
}


calculateAdjacency <- function(location.mat, window.mat){
  adjmat = Matrix::crossprod(location.mat, window.mat)
  Matrix::diag(adjmat) = 0  
  
  test = Matrix(sample(0:1,12,T), ncol=3, sparse=T)
  termfreq = Matrix::colSums(as(location.mat, 'dgCMatrix'))
  list(adjmat=adjmat, termfreq=termfreq)
}
