
#' Create a co-occurence network
#' 
#' @param x either a document-term matrix (DTM) or a tokenlist. 
#' @param measure the measure to calcualte adjacency. Currently supports cosine and conditional probability
#' @return A graph in the Igraph format in which edges represent the adjacency of terms
#' @export
coOccurenceNetwork <- function(x, measure='cooccurence', doc.col=getOption('doc.col','doc_id'), position.col=getOption('position.col','position'), word.col=getOption('word.col','word')){
  if('data.frame' %in% class(x)) {
    verifyTokenlistColumns(tokenlist, doc.col, position.col, word.col)
    x = tokenlistToDTM(x, doc.col, position.col, word.col)
  }
  if('DocumentTermMatrix' %in% class(x)) x = dtmToSparseMatrix(x)

  x = as(as(x, 'dgCMatrix'), 'dgTMatrix')
  if(measure == 'cosine') {
    mat = as(getCosine(x), 'dgCMatrix')
    g = graph.adjacency(mat, mode='upper', diag=F, weighted=T)
  }
  if(measure == 'cooccurence') {
    mat = getCoOccurence(x)
    g = graph.adjacency(mat, mode='upper', diag=F, weighted=T)
  }
  if(measure == 'conprob') {
    mat = getConditionalProbability(x)
    g = graph.adjacency(mat, mode='directed', diag=F, weighted=T)
  }
  V(g)$freq = col_sums(x)
  g = set.edge.attribute(g, measure, value=E(g)$weight)
  class(g) = c('semnet',class(g))
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

calculateAdjacency <- function(position.mat, window.mat){
  adj = Matrix::crossprod(position.mat, window.mat)
  g = graph.adjacency(adj, mode='directed', weighted=T, diag=F)
  V(g)$freq = col_sums(position.mat)
  E(g)$cooccurence = E(g)$weight
  E(g)$weight_pct = E(g)$weight / V(g)$freq[get.edgelist(g, names = F)[,1]]
  class(g) = c('semnet',class(g))
  g
}

aggCoOc <- function(x, position.mat, window.mat){
  cooc = position.mat[,x] & window.mat
  cooc = as(cooc, 'lgTMatrix')
  cooc = data.frame(x=x, y=cooc@j+1, context=cooc@i+1, weight=cooc@x)
  cooc = cooc[!cooc$x == cooc$y,]
  ddply(cooc, .(x,y,context), summarize, weight=sum(weight))
}

calculateAdjacencyPerContext <- function(position.mat, window.mat) {
  adj = ldply(1:ncol(position.mat), function(x) aggCoOc(x, position.mat, window.mat))
  adj$context = rownames(position.mat)[adj$context]
  adj$x = as.factor(colnames(position.mat)[adj$x])
  adj$y = as.factor(colnames(position.mat)[adj$y])
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
chi2 <- function(a,b,c,d, yates_correction=rep(F, length(a)), autocorrect=F){
  n = a+b+c+d
  sums = cbind(c1 = a+c, c2 = b+d, r1 = a+b, r2 = c+d)

  if(autocorrect){
    ## apply Cochrans criteria: no expected values below 1 and less than 20% of cells empty (which means none in a 2x2 design)
    ## if these are violated, use the yates_correction
    ## http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2041889/
    e = cbind(sums[,'c1'] / n, sums[,'c2'] / n)
    e = cbind(e * sums[,'r1'], e * sums[,'r2'])
    c1 = rowSums(e < 1) > 0          # at least one expected value below 1
    c2 = rowSums(sums < 5) > 0       # at least one cell below 5
    yates_correction = ifelse(c1 | c2, T, F)
  }
  x = a*d - b*c
  x = ifelse(yates_correction, abs(x) - n/2, x)
  chi = n*x^2 / (sums[,'c1'] * sums[,'c2'] * sums[,'r1'] * sums[,'r2'])
  ifelse(is.na(chi), 0, chi)
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

smooth_prob <- function(succes, n, vocabulary_size=NULL, smoothing_parameter=1){
  if(is.null(vocabulary_size)) {
    (word_occ + smoothing_parameter) / (n + smoothing_parameter)
  } else{
    (word_occ + smoothing_parameter) / (n + (vocabulary_size*smoothing_parameter))
  }
}

edgeChi2 <- function(g, wordfreq.col='freq', wordcooc.col='cooccurence', vocabulary=NULL){
  freq = get.vertex.attribute(g, wordfreq.col)
  cooc = get.edge.attribute(g, wordcooc.col)
  e = get.edgelist(g, names = F)
  if(is.null(vocabulary)) vocabulary = structure(freq, names=V(g)$name)

  y_if_x= cooc
  y_ifnot_x = freq[e[,2]] - cooc
  
  x_total_cooc = tapply(cooc, e[,1], sum)
  n_if_x = x_total_cooc[e[,1]]
  n_ifnot_x = sum(x_total_cooc) - n_if_x
  
  E(g)$chi =  chi2(y_if_x, y_ifnot_x, 
                   n_if_x - y_if_x, n_ifnot_x - y_ifnot_x,
                   autocorrect = T)
  
  
  E(g)$chi.p = 1-pchisq(E(g)$chi, 1)
  g
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
    class(g) = c('semnet',class(g))
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
    class(g) = c('semnet',class(g))
    return(g)
  } else return(d[,c('x','y', 'n', 'odds_ratio','smooth_odds_ratio','chi','p')])
}

edgeWeightStatistics <- function(g, edge.coocc.var='weight', vertex.freq.var='freq'){
  
}
