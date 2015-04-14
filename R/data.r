#' Document-term matrix based on wikinews articles on Iraq
#' 
#' @docType data
#' @keywords datasets
#' @name wikinews_iraq
#' @usage data(wikinews_iraq)
#' @format iraq.dtm: A tm document-term matrix of wikinews articles, and iraq.meta: a data frame containing meta-information about the articles in the dtm
NULL

#' Tokens from Bush and Obama's State of the Union addresses
#' 
#' @docType data
#' @keywords datasets
#' @name sotu
#' @usage data(sotu)
#' @format sotu.tokens: A dataframe of tokens produced by using Stanford CoreNLP to lemmatize the State of the Union speeches; and sotu.meta the metadata for these documents
NULL


#' Document metadata from Bush and Obama's State of the Union addresses
#' 
#' @docType data
#' @keywords datasets
#' @name sotu
#' @usage data(sotu)
#' @format sotu-meta: A dataframe of document metadata matching the tokens from sotu.tokens
NULL


#' Document-term matrix containing the nouns, names and adjectives from Bush and Obama's State of the Union addresses
#' 
#' @docType data
#' @keywords datasets
#' @name sotu
#' @usage data(sotu)
#' @format sotu.dtm: A document-term matrix based on sotu.tokens, containing only names, nouns, and adjectives
NULL