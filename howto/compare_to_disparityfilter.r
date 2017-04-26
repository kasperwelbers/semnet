# devtools::install_github('kasperwelbers/semnet')
# devtools::install_github('alessandrobessi/disparityfilter')
library(semnet) 
library(disparityfilter) 

######################## undirected networks. Gives the same results
set.seed(1)
g <- sample_pa(n = 500, m = 5, directed = FALSE)
E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)

## same edges and alpha, though from/to can be switched (undirected network)
disparityfilter::backbone(g)
get.data.frame(semnet::getBackboneNetwork(g, delete.isolates = F)) 


######################## directed network. Gives the same results
set.seed(1)
g <- sample_pa(n = 500, m = 5, directed = TRUE)
E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)

disparityfilter::backbone(g) 
get.data.frame(semnet::getBackboneNetwork(g, delete.isolates = F))


######################### directed network with parallel edges. Does not work with semnet
set.seed(1)
g <- sample_pa(n = 500, m = 5, directed = TRUE, algorithm = 'psumtree-multiple')
E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)

disparityfilter::backbone(g) 
get.data.frame(semnet::getBackboneNetwork(g, delete.isolates = F))
