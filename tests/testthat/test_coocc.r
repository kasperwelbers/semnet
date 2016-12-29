#library(semnet)

test_that("Document cooccurence", {
  
  ## DTM BASED COOCCURENCE ##
  data(simple_dtm, package='semnet')
  
  ## coOccurence freq
  g = coOccurenceNetwork(dtm)
  expect_equal(get.data.frame(g)$from, c("nuclear", "nuclear", "nuclear", "waste", "weapons", "weapons", "bad", "energy"))
  ## coOccurence cosine
  g = coOccurenceNetwork(dtm, 'cosine')
  expect_equal(round(E(g)$weight, 3), round(c(0.4082483, 0.4082483, 0.4082483, 0.5000000, 0.5000000, 0.7071068, 0.7071068, 0.7071068),3))
})

test_that("Windowed cooccurence", {
  ## TOKENLIST BASED COOCCURENCE WITH SLIDING WINDOW ##
  data(sotu_tokenlist, package='semnet')
  sotu.tokens = sotu.tokens[sotu.tokens$pos1 %in% c('N','M','A'),]
  sotu.tokens = sotu.tokens[sample(1:nrow(sotu.tokens)),]
  
  ## windowed cooccurence (using non-default columnnames)
  g = windowedCoOccurenceNetwork(sotu.tokens,
                                 window.size=20,
                                 doc.col='aid',
                                 position.col = 'id',
                                 word.col = 'lemma')
  expect_equal(get.data.frame(g)$from[1:10], c("task","basic","bargain","country","idea","people","sure","business","work","Congress"))
  expect_equal(E(g)$weight[1:10], c(3,1,1,2,1,1,1,2,1,1))
  
  ## NETWORK ANALYSIS AND VISUALIZATION ##
  
  ## backbone extraction
  g_backbone = getBackboneNetwork(g, alpha=0.0001, max.vertices=100)
  expect_equal(ecount(g_backbone), 255)
  
  ## pretty plotting
  gs = plot_semnet(g_backbone)
  expect_equal(unique(V(gs)$color)[1:10], c("#FF6666", "#FFB666", "#FF8E66", "#D066FF", "#FFDE66", "#66FF9B", "#A9FF66", "#66C3FF", "#D0FF66", "#F8FF66"))
  expect_equal(round(V(gs)$size[1:10],3), c(11.456, 13.544, 15.000,  9.071,  8.303,  9.471,  7.194,  7.637, 12.072,  5.673))
})


