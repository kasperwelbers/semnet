test_that("Query search works", {
  library(semnet)
  library(testthat)
  text = c('Renewable fuel is better than fossil fuels!',
           'A fueled debate about fuel',
           'Mark_Rutte is simply Rutte')
  tokens = asTokenlist(text)
  
  tokens = tokens[sample(1:nrow(tokens)),]
  setTokenlistColnames()

  ## simple indicator only
  hits = searchQuery(tokens, indicator = 'fuel')
  expect_equal(as.character(hits$word), c('fuel','fuel'))

  ## two indicators
  hits = searchQuery(tokens, indicator = 'fuel fuels')
  expect_equal(as.character(hits$word), c('fuel','fuels','fuel'))

  ## indicator with wildcard
  hits = searchQuery(tokens, indicator = 'fuel*')
  expect_equal(as.character(hits$word), c('fuel','fuels','fueled','fuel'))

  ## indicator and condition
  hits = searchQuery(tokens, indicator = 'fuel*', condition = 'renewable green clean')
  expect_equal(as.character(hits$word), c('fuel','fuels'))

  ## condition with default.window
  hits = searchQuery(tokens, indicator = 'fuel*', condition = 'renewable green clean', default.window = 2)
  expect_equal(as.character(hits$word), c('fuel'))

  ## condition once parameter
  hits_f = searchQuery(tokens, indicator = 'rutte', condition = 'mark~2')
  hits_t = searchQuery(tokens, indicator = 'rutte', condition = 'mark~2', condition_once = T)
  expect_equal(as.character(hits_f$word), c('mark rutte'))
  expect_equal(as.character(hits_t$word), c('mark rutte','rutte'))

  ## multiple queries
  queries = data.frame(code=c('renewable fuel', 'mark rutte', 'debate'),
                       indicator=c('fuel*', 'rutte', 'debate'),
                       condition = c('renewable green clean', 'mark~2', ''))
  hits = searchQueries(tokens, queries, condition_once=c(F,T,F))
  expect_equal(as.character(hits$word), c('fuel','fuels','mark rutte', 'rutte', 'debate'))

  ## batchsize = 1
  hits = searchQueries(tokens, queries, condition_once=c(F,T,F), batchsize=1)
  expect_equal(as.character(hits$word), c('fuel','fuels','mark rutte', 'rutte', 'debate'))

  ## code queries
  code = codeQueries(tokens, queries, condition_once=c(F,T,F))
  expect_equal(as.numeric(table(code)), c(12,1,2,2))
})
