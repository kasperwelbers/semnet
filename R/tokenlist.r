#' Create a tokenlist data.frame
#'
#' A tokenlist is a data.frame in which rows represent the tokens of a text (e.g., words, lemma, ngrams). This function creates a tokenlist that is ordered by document ('doc_id' column) and the position of the token in the text ('position' column).
#'
#' The tokenization is taken care of by the tokenize function of the quanteda package. Additional arguments (...) are passed to the tokenize function.
#'
#' The default column names for the tokenlist are "doc_id", "position" and "word". Functions in semnet where the tokenlist should be given as an argument assume that these column names are used. 
#' If alternative columnnames are prefered, these can be specified in two ways. First, the defaults can be set when calling a function using the doc.col, position.col and word.col parameters. Second, defaults can be set globally by using the setTokenlistColnames() function.
#'
#' @param x An object that can be transformed into a tokenlist object. This can be 1) a list of the tokenizedTexts class (quanteda). 2) A data.frame with document_id, position and word columns (see above for explanation of columnnames). Or 3) a character vector, in which case the tokenize function of the quanteda package is used.
#' @param doc_id If the input is a tokenizedTexts list or character vector, the doc_id vector can be given to define document ids (otherwise, the list or vector indices are used)
#' @param doc.col The name of the document_id column. Defaults to "doc_id", unless a global default is specified using setTokenlistColnames()
#' @param position.col The name of the column giving the position in a document. Defaults to "position", unless a global default is specified using setTokenlistColnames()
#' @param word.col The name of the column containing the token text. Defaults to "word", unless a global default is specified using setTokenlistColnames()
#' @param ... If x is a character vector, additional arguments will be passed to the tokenize function of the quanteda package
#'
#' @return
#' @export
asTokenlist <- function(x, doc_id=NULL, language='english', use_stemming=F, lowercase=T, remove_stopwords=F, ngrams=1, doc.col = getOption('doc.col','doc_id'), position.col = getOption('position.col','position'), word.col=getOption('word.col','word'), ...){
  if(is(x, 'factor')) x = as.character(x)
  if(is(x, 'character')) x = quanteda::tokenize(x, ngrams=ngrams, ...)
  if(is(x, 'tokenizedTexts')) x = tokenizedTextsToDf(x, doc_id, doc.col, position.col, word.col)

  ## at this point, x should be a data.frame with columnnames matching doc.col, position.col and word.col
  if(!is(x, 'data.frame')) stop('Tokenlist argument should be one of the following:\n- tokenizedTexts class (quanteda)\n- corpus class (quanteda)\n- tokenlist class (semnet)\n- a data.frame with document_id, position and word columns (see the *.col arguments)\n- a character vector')
  if(!class(x[,word.col]) == 'factor') x[,word.col] = as.factor(x[,word.col])
  verifyTokenlistColumns(x, doc.col, position.col, word.col)
  
  ## pre-processing
  if(lowercase) levels(x[,word.col]) = tolower(levels(x[,word.col]))
  if(remove_stopwords) x = x[!x[,word.col] %in% stopwords(language),]
  if(use_stemming) levels(x[,word.col]) = quanteda::wordstem(levels(x[,word.col]), language=language)
  
  class(x) = c('tokenlist',class(x))
  x
}

tokenizedTextsToDf <- function(x, doc_id, doc.col, position.col, word.col){
  if(is.null(doc_id)) doc_id = 1:length(x)
  doclen = sapply(x, length)
  
  filter = doclen > 0 ## ignore document with length 0
  x = data.frame(doc_id = rep(doc_id[filter], doclen[filter]),
                 position = unlist(sapply(doclen[filter], function(x) 1:x, simplify = F)),
                 word = unlist(x[filter]))
  colnames(x) = c(doc.col, position.col, word.col)
  x
}

verifyTokenlistColumns <- function(x, doc.col = getOption('doc.col','doc_id'), position.col = getOption('position.col','position'), word.col=getOption('word.col','word')){
  default = c(doc.col, position.col, word.col)
  wrong = default[!default %in% colnames(x)]
  
  if(length(wrong) > 0){
    wrong = paste(sprintf('"%s"', wrong), collapse=', ')
    mes = sprintf('The following specified or default columnnames do not occur in the tokenlist: %s.\nPlease view ?asTokenlist() for instructions.', wrong)
    stop(mes)
  }
}

#' Transform a tokenlist into a document-term matrix (DTM)
#'
#' @param tokenlist 
#' @param doc.col 
#' @param position.col 
#' @param word.col 
#'
#' @return a DTM, as a sparse matrix in the dgTMatrix class 
#' @export
tokenlistToDTM <- function(tokenlist, doc.col = getOption('doc.col','doc_id'), position.col = getOption('position.col','position'), word.col=getOption('word.col','word')){
  udoc = unique(tokenslist[,doc.col])
  uterm = unique(tokens[,word.col])
  dtm = spMatrix(length(udoc), length(uterm), match(tokenlist[,doc.col], udoc), match(tokenlist[,word.col], uterm), rep(1, nrow(tokenlist)))
  dimnames(dtm) = list(udoc, uterm)
  as(as(dtm,'dgCMatrix'), 'dgTMatrix')
}

