##### windowed adjacency functions #####

#' A sliding window approach to calculate the co-occurence of words
#' 
#' @param tokenlist 
#' @param window.size The distance within which words should occur from each other to be counted as a co-occurence.
#' @param output.per.context Logical. If True, co-occurences are reported per context (beware that this takes longer and can lead to huge output)
#' @param direction a string indicating whether only the left ('<') or right ('>') side of the window, or both ('<>'), should be used. 
#'
#' @return An edgelist (data.frame) with columns x, y and weight, in which weight represents the number of times y occured within a [window.size] word distance from x. If output.per.context is True, co-occurences are reported per context, and the edgelist has an additional context column.
#' @export
windowedCoOccurenceNetwork <- function(tokenlist, window.size=10, output.per.context=F, direction='<>', doc.col=getOption('doc.col','doc_id'), position.col=getOption('position.col','position'), word.col=getOption('word.col','word')){
  verifyTokenlistColumns(tokenlist, doc.col, position.col, word.col)
  
  if(min(tokenlist[,position.col]) == 0) tokenlist[,position.col] = tokenlist[,position.col] + 1 # if indexing starts at 0, set to 1
  mat = wordWindowOccurence(tokenlist, window.size, direction, doc.col=doc.col, position.col=position.col, word.col=word.col)
  if(output.per.context) {
    calculateAdjacencyPerContext(mat$position.mat, mat$window.mat)
  } else {
    calculateAdjacency(mat$position.mat, mat$window.mat)
  }
}

getWindowMatrix <- function(position, context, term, window.size, return_i_filter, presorted){
  position = globalPosition(position, context, window.size=window.size, presorted)
  shifts = -window.size:window.size
  
  if(!is.null(return_i_filter)) return_i_filter = position[return_i_filter] ## make sure the return_i_filter uses the transformed indices
  
  terms = unique(term)
  term_index = match(term, terms)
  m = positionMatrix(i=position, j=term_index, shifts, distance.as.value=T, return_i_filter=return_i_filter)
  colnames(m) = terms
  m
}

localPosition <- function(position, context, presorted=F){
  if(!presorted){
    ord = order(context, position)
    position = position[ord]
    context = context[ord]
  }
  newcontext = which(!duplicated(context))
  repeat_multiplier = c(newcontext[-1], length(context)+1) - newcontext
  context_start = rep(position[newcontext], repeat_multiplier)
  position = (position - context_start) + 1
  if(!presorted) position = position[match(1:length(position), ord)] 
  position
}

globalPosition <- function(position, context, window.size=NA, presorted=F){
  ## makes the word position counter global with dummy positions between contexts to prevent overlapping windows (so it can be used as an index).
  ## this way, overlapping word windows can be calculated for multiple documents within a single matrix.
  ## position and context need to be sorted on order(context,position). If this is already the case, presorted=T can be used for a speed up
  if(!presorted){
    ord = order(context, position)
    position = position[ord]
    context = context[ord]
  }
  
  ## first, make sure position is local and starts at 1 for each context (otherwise things get very slow)
  position = localPosition(position, context, presorted=T)
  
  if(min(position) == 0) position = position + 1 ## position will be treated as an index, so it cannot be zero in r where an index starts at 1 (and some parsers start indexing at zero)
  
  if(!length(unique(context)) == 1) {
    newcontext = which(!duplicated(context)) # where does a new context start
    
    context.max = position[newcontext-1] # the highest value of each context
    if(!is.na(window.size)) context.max = context.max + window.size # increase the highest value of each context with window.size to make sure windows of different contexts do not overlap.
    multiplier_scores = cumsum(c(0,context.max)) # the amount that should be added to the position at the start of each context
    
    repeat_multiplier = c(newcontext[-1], length(position)+1) - newcontext # the number of times the multiplier scores need to be repeated to match the position vector
    multiplier_vector = rep(multiplier_scores, repeat_multiplier)
    position = position + multiplier_vector
  }
  if(!presorted) position = position[match(1:length(position), ord)] 
  position
}

positionMatrix <- function(i, j, shifts=0, count.once=T, distance.as.value=F, return_i_filter=NULL){
  mat = spMatrix(max(i), max(j))
  
  shifts = shifts[order(abs(shifts))] # order from 0 to higher (required if distance.as.value = T)
  for(shift in shifts){
    i_shift = i + shift
    
    if(!is.null(return_i_filter)) {
      select = i_shift %in% return_i_filter
    } else {
      select = i_shift > 0 & i_shift <= max(i)
    }
    
    if(sum(select) == 0) next
    
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
  mat = mat[i,,drop=F]
  mat = as(mat, 'dgCMatrix')
  if(count.once) mat@x[mat@x>0] = 1
  mat
}

#' Creates a matrix where rows match the tokens, columns represent the queries in query_regex, and values represent the word distance to each query hit
#'
#' @param tokens
#' @param query_regex
#' @param word.col
#' @param presorted
#' @param default.window
#'
#' @return a sparse matrix
#' @export
getQueryMatrix <- function(tokens, query_regex, doc.col=getOption('doc.col','doc_id'), position.col=getOption('position.col','position'), word.col=getOption('word.col','word'), presorted=F, default.window=NA, return_i=NULL){
  if(!'window' %in% colnames(query_regex)) query_regex$window = default.window
  
  ## replace terms that do not occur in the queries with a dummy (NA). Skip if there are already NA's in the data
  if(sum(is.na(tokens[,word.col])) == 0) {
    tokens[,word.col] = ifelse(tokenGrepl(query_regex$regex, tokens[,word.col]), as.character(tokens[,word.col]), NA)
  }
  
  ## get query matrix; separately for query terms with a given word distance and query terms at the article level
  document_level = is.na(query_regex$window) | query_regex$window == 'd'
  if(mean(document_level) == 1) qm = documentOccurenceQueryMatrix(tokens, query_regex[document_level,], doc.col, position.col, word.col)
  if(mean(document_level) == 0) qm = wordDistanceQueryMatrix(tokens, query_regex[!document_level,], doc.col, position.col, word.col, return_i, presorted)
  if(!mean(document_level) %in% c(0,1)) qm = cbind(wordDistanceQueryMatrix(tokens, query_regex[!document_level,], doc.col, position.col, word.col, return_i, presorted),
                                                   documentOccurenceQueryMatrix(tokens, query_regex[document_level,], doc.col, position.col, word.col))
  
  if(!is.null(return_i)) qm = qm[return_i,,drop=F]
  qm
}

wordDistanceQueryMatrix <- function(tokens, query_regex, doc.col, position.col, word.col, return_i_filter=NULL, presorted=F){
  if(nrow(query_regex) == 0) return(NULL)
  query_regex$window = as.numeric(query_regex$window) + 1 # Plus 1, because in the window matrix 1 indicates no distance (because zero already indicates no presence [because this keeps the matrix sparse])
  
  m = getWindowMatrix(position=tokens[,position.col],
                      context=tokens[,doc.col],
                      term = tokens[,word.col],
                      window.size = max(query_regex$window),
                      return_i_filter = return_i_filter,
                      presorted=presorted)
  
  ## create the rows and columns for the query matrix by looking in which rows one of the terms that matches the regex is TRUE.
  getWindowQueryHits <- function(j, query_regex, m){
    query_m = m[,tokenGrepl(query_regex$regex[j], colnames(m)),drop=F]
    hits = Matrix::rowSums(query_m > 0 & query_m <= query_regex$window[j]) > 0 # at least one column should have a value between 1 and its window size
    if(sum(hits) == 0) return(NULL)
    data.frame(i = which(hits), j = j)
  }
  
  
  qm = ldply(1:nrow(query_regex), getWindowQueryHits, query_regex=query_regex, m=m)
  qm = spMatrix(nrow(tokens), nrow(query_regex), qm$i, qm$j, rep(T, nrow(qm)))
  colnames(qm) = query_regex$term
  qm
}

documentOccurenceQueryMatrix <- function(tokens, query_regex, doc.col, position.col, word.col){
  if(nrow(query_regex) == 0) return(NULL)
  
  udoc = unique(tokens[,doc.col])
  uterm = unique(tokens[,word.col])
  doc_i = match(tokens[,doc.col], udoc)
  term_i = match(tokens[,word.col], uterm)
  m = spMatrix(length(udoc), length(uterm), doc_i, term_i, rep(1,nrow(tokens)))
  colnames(m) = uterm
  
  ## create the rows and columns for the query matrix by looking in which rows one of the terms that matches the regex is TRUE.
  getQueryHits <- function(j, query_regex, m){
    query_m = m[,tokenGrepl(query_regex$regex[j], colnames(m)),drop=F]
    hits = Matrix::rowSums(query_m) > 0 # matrix value has to be higher than 0, but not higher than the window size
    if(sum(hits) == 0) return(NULL)
    data.frame(i = which(hits), j = j)
  }
  
  qm = ldply(1:nrow(query_regex), getQueryHits, query_regex=query_regex, m=m)
  
  qm = spMatrix(nrow(tokens), nrow(query_regex), qm$i, qm$j, rep(T, nrow(qm)))
  colnames(qm) = query_regex$term
  
  qm[doc_i,,drop=F] ## repeat rows (i.e. documents) to match the token list input
}


#' Gives the window in which a term occured in a matrix.
#' 
#' This function returns the occurence of words (position.matrix) and the window of occurence (window.matrix). This format enables the co-occurence of words within sliding windows (i.e. word distance) to be calculated by multiplying position.matrix with window.matrix. 
#' 
#' @param position An integer vector giving the position of terms in a given context (e.g., document, paragraph, sentence) 
#' @param term A character vector giving the terms
#' @param context A vector giving the context in which terms occur (e.g., document, paragraph, sentence)
#' @param window.size The distance within which words should occur from each other to be counted as a co-occurence.
#' @param direction a string indicating whether only the left ('<') or right ('>') side of the window, or both ('<>'), should be used. 
#' @return A list with two matrices. position.mat gives the specific position of a term, and window.mat gives the window in which each word occured. The rows represent the position of a term, and matches the input of this function (position, term and context). The columns represents terms.
#' @export
wordWindowOccurence <- function(tokenlist, window.size=3, direction='<>', distance.as.value=F, doc.col=getOption('doc.col','doc_id'), position.col=getOption('position.col','position'), word.col=getOption('word.col','word')){
  tokenlist = tokenlist[!is.na(tokenlist[,word.col]),]
  
  ## order on context and position
  ord = order(tokenlist[,doc.col], tokenlist[,position.col])
  tokenlist = tokenlist[ord,]
  
  tokenlist[,position.col] = globalPosition(tokenlist[,position.col], tokenlist[,doc.col], window.size=window.size, presorted = T)

  if(direction == '<') shifts = -window.size:0
  if(direction == '<>') shifts = -window.size:window.size
  if(direction == '>') shifts = 0:window.size
  
  terms = unique(tokenlist[,word.col])
  term_index = match(tokenlist[,word.col], terms)
  
  position.mat = positionMatrix(tokenlist[,position.col], term_index, 0)
  window.mat = positionMatrix(tokenlist[,position.col], term_index, shifts, distance.as.value=distance.as.value)
  
  colnames(position.mat) = colnames(window.mat) = terms
  rownames(position.mat) = rownames(window.mat) = tokenlist[!duplicated(tokenlist[,position.col]), doc.col]

  ## return original order
  if(!identical(ord, 1:nrow(position.mat))){
    inverse_ord = match(1:nrow(position.mat), ord)
    position.mat = position.mat[inverse_ord,] 
    window.mat = window.mat[inverse_ord,] 
  }
  
  list(position.mat=position.mat, window.mat=window.mat)
}