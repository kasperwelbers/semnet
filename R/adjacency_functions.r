

#' Create an adjacency graph from a document term matrix
#' 
#' @param dtm a document term matrix, either or not in the dtm format from the `tm` package
#' @param measure the measure to calcualte adjacency. Currently supports cosine and conditional probability
#' @return A graph in the Igraph format in which edges represent the adjacency of terms
#' @export
coOccurenceNetwork <- function(dtm, measure='cooccurence'){
  if('DocumentTermMatrix' %in% class(dtm)) dtm = dtmToSparseMatrix(dtm)
  dtm = as(dtm, 'dgCMatrix')
  if(measure == 'cosine') {
    mat = getCosine(dtm)
    Matrix::diag(mat) = 0 # ignore loops
    g = graph.adjacency(mat, mode='upper', diag=F, weighted=T)
  }
  if(measure == 'cooccurence') {
    mat = getCoOccurence(dtm)
    g = graph.adjacency(mat, mode='upper', diag=F, weighted=T)
  }
  if(measure == 'conprob') {
    mat = getConditionalProbability(dtm)
    g = graph.adjacency(mat, mode='directed', diag=F, weighted=T)
  }
  V(g)$freq = col_sums(dtm)
  #edgelist = Matrix::which(mat>0, arr.ind=T)
  #edgelist = data.frame(x=rownames(mat)[edgelist[,1]], y=colnames(mat)[edgelist[,2]], weight=mat@x[mat@x>0])
  #edgelist
  g
}

dtmToSparseMatrix <- function(dtm){
  sm = spMatrix(nrow(dtm), ncol(dtm), dtm$i, dtm$j, dtm$v)
  rownames(sm) = rownames(dtm)
  colnames(sm) = colnames(dtm)
  sm
}

getCoOccurence <- function(m1, m2=m1){
  m1@x[m1@x > 0] = 1
  m2@x[m2@x > 0] = 1
  mat = Matrix::crossprod(m1,m2) 
  mat[is.na(mat)] = 0
  mat
}

getConditionalProbability <- function(m1, m2=m1){
  m1@x[m1@x > 0] = 1
  m2@x[m2@x > 0] = 1
  mat = Matrix::crossprod(m1,m2)/colSums(m1) 
  mat[is.na(mat)] = 0
  mat
}

getCosine <- function(m1, m2=m1){
  norm.x = sqrt(colSums(m1^2))
  norm.y = sqrt(colSums(m2^2))
  mat = Matrix::crossprod(m1,m2)
  mat = mat / Matrix::tcrossprod(norm.x, norm.y)
  mat[is.na(mat)] = 0
  mat
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
  mat = as(mat, 'dgCMatrix')
  if(count.once) mat@x[mat@x>0] = 1
  mat
}

#' Gives the window in which a term occured in a matrix.
#' 
#' This function returns the occurence of words (location.matrix) and the window of occurence (window.matrix). This format enables the co-occurence of words within sliding windows (i.e. word distance) to be calculated by multiplying location.matrix with window.matrix. 
#' 
#' @param location An integer vector giving the position of terms in a given context (e.g., document, paragraph, sentence) 
#' @param term A character vector giving the terms
#' @param context A vector giving the context in which terms occur (e.g., document, paragraph, sentence)
#' @param window.size The distance within which words should occur from each other to be counted as a co-occurence.
#' @param two.sided Logical. If false, it is only counted how often a word occured `after` another word within the given window size
#' @return A list with two matrices. location.mat gives the specific location of a term, and window.mat gives the window in which each word occured. The rows represent the location of a term, and matches the input of this function (location, term and context). The columns represents terms.
#' @export
wordWindowOccurence <- function(location, term, context, window.size=3, two.sided=T){
  nas = is.na(term)
  if (any(nas)) {
    term = term[!nas]
    location = location[!nas]
    context = context[!nas]
  }
  
  ord = order(context, location)
  location = location[ord]
  term = term[ord]
  context = context[ord]
  
  location = stretchLocation(location,context,window.size=window.size)
  shifts = if(two.sided) -window.size:window.size else 0:window.size
  terms = unique(term)
  term_index = match(term, terms)
  
  location.mat = locationMatrix(location, term_index, 0)
  window.mat = locationMatrix(location, term_index, shifts)
  
  colnames(location.mat) = colnames(window.mat) = terms
  rownames(location.mat) = rownames(window.mat) = context

  
  list(location.mat=location.mat, window.mat=window.mat)
  
  #calculateAdjacency(location.mat, window.mat)
}

#' A sliding window approach to calculate the co-occurence of words
#' 
#' @param location An integer vector giving the position of terms in a given context (e.g., document, paragraph, sentence) 
#' @param term A character vector giving the terms
#' @param context A vector giving the context in which terms occur (e.g., document, paragraph, sentence)
#' @param window.size The distance within which words should occur from each other to be counted as a co-occurence.
#' @param output.per.context Logical. If True, co-occurences are reported per context (beware that this takes longer and can lead to huge output)
#' @param two.sided Logical. If false, it is only counted how often a word occured `after` another word within the given window size
#' @return An edgelist (data.frame) with columns x, y and weight, in which weight represents the number of times y occured within a [window.size] word distance from x. If output.per.context is True, co-occurences are reported per context, and the edgelist has an additional context column.
#' @export
windowedCoOccurenceNetwork <- function(location, term, context, window.size=10, output.per.context=F, two.sided=T){
  if(min(location) == 0) location = location + 1 # if indexing starts at 0, set to 1
  mat = wordWindowOccurence(location, term, context, window.size, two.sided)
  if(output.per.context) {
    calculateAdjacencyPerContext(mat$location.mat, mat$window.mat)
  } else {
    calculateAdjacency(mat$location.mat, mat$window.mat)
  }
}

calculateAdjacency <- function(location.mat, window.mat){
  adj = Matrix::crossprod(location.mat, window.mat)
  g = graph.adjacency(adj, mode='directed', weighted=T, diag=F)
  V(g)$freq = col_sums(location.mat)
  g
  #edgelist = Matrix::which(adj>0, arr.ind=T)
  #edgelist = data.frame(x=rownames(adj)[edgelist[,1]], y=colnames(adj)[edgelist[,2]], weight=adj@x[adj@x>0])
  #edgelist
}

aggCoOc <- function(x, location.mat, window.mat){
  cooc = location.mat[,x] & window.mat
  cooc = as(cooc, 'lgTMatrix')
  cooc = data.frame(x=x, y=cooc@j+1, context=cooc@i+1, weight=cooc@x)
  cooc = cooc[!cooc$x == cooc$y,]
  ddply(cooc, .(x,y,context), summarize, weight=sum(weight))
}

calculateAdjacencyPerContext <- function(location.mat, window.mat) {
  adj = ldply(1:ncol(location.mat), function(x) aggCoOc(x, location.mat, window.mat))
  adj$context = rownames(location.mat)[adj$context]
  adj$x = as.factor(colnames(location.mat)[adj$x])
  adj$y = as.factor(colnames(location.mat)[adj$y])
  adj
}

#' Merge a list of matrices
#' 
#' Merge matrices, by standradizing column and row dimensions and using the Reduce function to merge them.
#' 
#' @param matrices A list of matrices
#' @param reduce_func The function to pass to Reduce
#' @return A single matrix
#' @export
mergeMatrices <- function(matrices, reduce_func='+'){
  rowdim = unique(llply(matrices, function(x) rownames(x)))[[1]]
  coldim = unique(llply(matrices, function(x) colnames(x)))[[1]]
  matrices = llply(matrices, function(x) setMatrixDims(x, rowdim, coldim))
  mat = Reduce(reduce_func, matrices)
  rownames(mat) = rowdim
  colnames(mat) = coldim
  mat
}

setMatrixDims <- function(mat, rowdim, coldim){
  ## set rows and columns of matrix to given rowdim and coldim vectors. (usefull for making different matrices identical so that they can be summed up)  
  mat = as(mat, 'dgTMatrix')
  d = data.frame(i=rownames(mat)[mat@i+1], j=colnames(mat)[mat@j+1], v=mat@x) 
  d = d[d$i %in% rowdim & d$j %in% coldim,]
  d$i = match(d$i, rowdim)
  d$j = match(d$j, coldim)
  spMatrix(length(rowdim), length(coldim), d$i, d$j, d$v)
}

