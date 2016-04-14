
#' Create an adjacency graph from a document term matrix
#' 
#' @param dtm a document term matrix, either or not in the dtm format from the `tm` package
#' @param measure the measure to calcualte adjacency. Currently supports cosine and conditional probability
#' @return A graph in the Igraph format in which edges represent the adjacency of terms
#' @export
coOccurenceNetwork <- function(dtm, measure='cooccurence'){
  if('DocumentTermMatrix' %in% class(dtm)) dtm = dtmToSparseMatrix(dtm)
  dtm = as(as(dtm, 'dgCMatrix'), 'dgTMatrix')
  if(measure == 'cosine') {
    mat = as(getCosine(dtm), 'dgCMatrix')
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
  mat = Matrix::crossprod(m1,m2)/Matrix::colSums(m1) 
  mat[is.na(mat)] = 0
  mat
}

getCosine <- function(m1, m2=NULL){
  norm = sqrt(Matrix::colSums(m1^2))
  m1@x = m1@x / norm[m1@j+1]  
  if(!is.null(m2)){
    norm = sqrt(Matrix::colSums(m2^2))
    m2@x = m2@x / norm[m2@j+1]
    cp = Matrix::crossprod(m1,m2) 
  } else cp = Matrix::crossprod(m1)
  cp
}

##### windowed adjacency functions #####

stretchLocation <- function(location, context, window.size){
  ## makes the word location counter global, and adds dummy locations between contexts to prevent overlapping windows.
  ## this way, overlapping word windows can be calculated for multiple documents within a single matrix.
  ## location and context need to be sorted on order(context,location)!!
  newcontext = which(!duplicated(context)) # where does a new context start
  context.max = location[newcontext-1] + (window.size*2)
  multiplier_scores = cumsum(c(0,context.max)) # the amount that should be added to the location at the start of each context 
  repeat_multiplier = c(newcontext[2:length(newcontext)], length(location)+1) - newcontext # the number of times the multiplier scores need to be repeated to match the location vector
  multiplier_vector = rep(multiplier_scores, repeat_multiplier)
  return(location + multiplier_vector)
}

locationMatrix <- function(i, j, shifts=0, count.once=T, distance.as.value=F){
  mat = spMatrix(max(i), max(j))
  
  shifts = shifts[order(abs(shifts))] # order from 0 to higher (required if distance.as.value = T)
  for(shift in shifts){
    i_shift = i + shift
    select = i_shift > 0 & i_shift <= max(i)
    if(distance.as.value){
      mat = mat + spMatrix(nrow=max(i), ncol=max(j), i=i_shift[select], j=j[select], rep(abs(shift)+1, sum(select))) 
    } else{
      mat = mat + spMatrix(nrow=max(i), ncol=max(j), i=i_shift[select], j=j[select], rep(1, sum(select)))   
    }
  }
  
  if(distance.as.value){
    ## remove duplicates. since the stacked triples are ordered by shifts, this leaves the shortest distance to a term in case of duplicate cells
    count.once = F
    select = !duplicated(data.frame(mat@i, mat@j))
    mat = spMatrix(nrow(mat), ncol(mat), mat@i[select]+1, mat@j[select]+1, mat@x[select])
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
wordWindowOccurence <- function(location, term, context, window.size=3, two.sided=T, distance.as.value=F){
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
  window.mat = locationMatrix(location, term_index, shifts, distance.as.value=distance.as.value)
  
  colnames(location.mat) = colnames(window.mat) = terms
  rownames(location.mat) = rownames(window.mat) = context

  
  list(location.mat=location.mat, window.mat=window.mat)
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


#### using chi-square (experimental)

#' Compute the chi^2 statistic for a 2x2 crosstab containing the values
#' [[a, b], [c, d]]
chi2 <- function(a,b,c,d) {
  ooe <- function(o, e) {(o-e)*(o-e) / e}
  tot = 0.0 + a+b+c+d
  a = as.numeric(a)
  b = as.numeric(b)
  c = as.numeric(c)
  d = as.numeric(d)
  (ooe(a, (a+c)*(a+b)/tot)
   +  ooe(b, (b+d)*(a+b)/tot)
   +  ooe(c, (a+c)*(c+d)/tot)
   +  ooe(d, (d+b)*(c+d)/tot))
}

binom.coef.log <- function(n, k) {
  bcl_func <- function(n, k) sum(log(((n-k+1):n) / ((k-k+1):k)))
  mapply(bcl_func, n, k)
}
fishers.exact <- function(a,b,c,d){  
  n = a+b+c+d
  log_p = binom.coef.log(a+b,a) + binom.coef.log(c+d,c) - binom.coef.log(n,a+c)
  exp(log_p)
}

getContextOccurence <- function(m1, m2=m1){
  m2@x[Matrix::which(m2@x > 0)] = 1
  mat = Matrix::crossprod(m1,m2)
  mat[is.na(mat)] = 0
  mat
}

smooth_prob <- function(succes, n, vocabulary_siz=NULL, smoothing_parameter=1){
  if(is.null(vocabulary_size)) {
    (word_occ + smoothing_parameter) / (n + smoothing_parameter)
  } else{
    (word_occ + smoothing_parameter) / (n + (vocabulary_size*smoothing_parameter))
  }
}

#' Get wordassociations based on ratio and chi-2
#' 
#' Calculate wordassociations (x -> y) as the ratio of the p
#' 
#' @param dtm the main document-term matrix
#' @return A dataframe representing and edge.list
#' @export
chi2_wordassociations <- function(dtm, odds.ratio.thres=1, thres=0.05, measure='fisher', smoothing_parameter=1, use.wordcount=T, return.graph=T) {
  d = coOccurenceDistributionTable(dtm, use.wordcount)
  
  d$prob_y_if_x = (d$y_if_x / d$n_if_x)
  d$prob_y_ifnot_x = (d$y_ifnot_x / d$n_ifnot_x)
  d$prob_ratio = d$prob_y_if_x / d$prob_y_ifnot_x
  d$odds_ratio = (d$prob_y_if_x / (1 + d$prob_y_if_x)) / (d$prob_y_ifnot_x / (1 + d$prob_y_ifnot_x))
  d = d[d$odds_ratio > odds.ratio.thres,]

  
  d$chi = chi2(d$y_if_x, 
               d$y_ifnot_x, 
               d$n_if_x-d$y_if_x, 
               d$n_ifnot_x-d$y_ifnot_x)
  d$chi.p = 1-pchisq(d$chi, 1)
  
  d$fisher.p = fishers.exact(d$y_if_x, 
                             d$y_ifnot_x, 
                             d$n_if_x-d$y_if_x, 
                             d$n_ifnot_x-d$y_ifnot_x)
  
  d = d[d$chi.p < chi.p.thres,]  
  d = d[d$fisher.p < fisher.p.thres,]  
  
  d = d[!d$x == d$y,]
  d$x = colnames(dtm)[d$x]
  d$y = colnames(dtm)[d$y]
  d$n = d$y_if_x
  
  if(return.graph){

    g = graph.data.frame(d[,c('x','y')])
    E(g)$odds_ratio = d$odds_ratio
    E(g)$smooth_odds_ratio = d$odds_ratio
    E(g)$weight = d$odds_ratio / (1 + d$odds_ratio)
    
    V(g)$freq = col_sums(dtm)[match(V(g)$name, colnames(dtm))]
    return(g)
  } else return(d[,c('x','y', 'n', 'prob_ratio', 'odds_ratio','chi','chi.p', 'fisher.p')])
}

coOccurenceDistributionTable <- function(dtm, use.wordcount=T, smooth_param=1){
  if('DocumentTermMatrix' %in% class(dtm)) dtm = dtmToSparseMatrix(dtm)
  dtm = as(dtm, 'dgCMatrix')
  
  conocc = as(getContextOccurence(dtm), 'dgTMatrix')    
  d = data.frame(x=conocc@j+1, y=conocc@i+1) 
  d$y_if_x= conocc@x # how often did y occur in documents in which x occured
  d$y_ifnot_x = Matrix::colSums(dtm)[d$y] - d$y_if_x # how often did y occur in documents in which x did not occur
  d$n_if_x = Matrix::colSums(conocc)[d$x] # how many words occured in documents where x occured
  d$n_ifnot_x = sum(dtm) - d$n_if_x # how many words occured in documents where x did not occur
  
  d$y_if_x = d$y_if_x + smooth_param
  d$y_ifnot_x = d$y_ifnot_x + smooth_param
  
  if(use.wordcount){
    vocabulary_if_x = Matrix::colSums(conocc > 0)[d$x]
    vocabulary = Matrix::colSums(dtm) 
    vocabulary_ifnot_x = sapply(1:ncol(dtm), function(x) sum(vocabulary > Matrix::colSums(dtm[dtm[,x] > 0, ,drop=F])))
    vocabulary_ifnot_x = vocabulary_ifnot_x[d$x]
    d$n_if_x = d$n_if_x + vocabulary_if_x*smooth_param
    #d$n_ifnot_x = d$n_ifnot_x + ncol(dtm)*smooth_param
    d$n_ifnot_x = d$n_ifnot_x + vocabulary_ifnot_x*smooth_param
    
  } else {
    d$n_if_x = d$n_if_x + smooth_param
    d$n_ifnot_x = d$n_ifnot_x + smooth_param
  }
  d
}

#' Get wordassociations based on ratio and chi-2
#' 
#' Calculate wordassociations (x -> y) as the ratio of the p
#' 
#' @param dtm the main document-term matrix
#' @return A dataframe representing and edge.list
#' @export
chi2_wordassociations_alternative <- function(dtm, ratio.thres=1, p.thres=0.05, smoothing_parameter=1, return.graph=T) {
  if('DocumentTermMatrix' %in% class(dtm)) dtm = dtmToSparseMatrix(dtm)
  dtm = as(dtm, 'dgCMatrix')
  conocc = as(getContextOccurence(dtm), 'dgTMatrix')  
  
  d = data.frame(x=conocc@j+1, y=conocc@i+1) 
  d$y_if_x= conocc@x # how often did y occur in documents in which x occured
  d$y_ifnot_x = diag(conocc)[d$y] - d$y_if_x # how often did y occur in documents in which x did not occur
  d$n_if_x = colSums(conocc)[d$x] # how many words occured in documents where x occured
  d$n_ifnot_x = sum(dtm) - d$n_if_x # how many words occured in documents where x did not occur
  
  d$chi = chi2(d$y_if_x, 
               d$y_ifnot_x, 
               d$n_if_x-d$y_if_x, 
               d$n_ifnot_x-d$y_ifnot_x)
  d$p = 1-pchisq(d$chi, 1)
  d = d[d$p < p.thres,]  
  
  
  d$odds_ratio = (d$y_if_x/d$n_if_x) / (d$y_ifnot_x/d$n_ifnot_x)
  d = d[d$odds_ratio > ratio.thres,]
  d$smooth_odds_ratio = smooth_rel(d$y_if_x, d$n_if_x, nrow(conocc), smoothing_parameter) / smooth_rel(d$y_ifnot_x, d$n_ifnot_x, nrow(conocc), smoothing_parameter)
  
  d = d[!d$x == d$y,]
  d$x = colnames(conocc)[d$x]
  d$y = colnames(conocc)[d$y]
  d$n = d$y_if_x
  
  if(return.graph){
    g = graph.data.frame(d[,c('x','y')])
    E(g)$odds_ratio = d$odds_ratio
    E(g)$smooth_odds_ratio = d$smooth_odds_ratio
    E(g)$weight = d$smooth_odds_ratio / (1 + d$smooth_odds_ratio)
    
    V(g)$freq = colSums(dtm)[match(V(g)$name, colnames(dtm))]
    return(g)
  } else return(d[,c('x','y', 'n', 'odds_ratio','smooth_odds_ratio','chi','p')])
}
