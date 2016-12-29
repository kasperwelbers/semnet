#cbfile = read.csv('~/Dropbox/syntax/data/cb_jan.csv')
codebookFromColumnTree <- function(cb, code_columns=grep('[0-9]+', colnames(cb), value=T)){
  code_columns = cb[,code_columns,drop=F]
  code_columns[is.na(code_columns)] = ''
  filter = !apply(code_columns, 1, paste, collapse='') == ''
  
  cb = data.frame(code = apply(code_columns[filter,], 1, paste, collapse=''),
                  level = apply(code_columns[filter,], 1, function(x) min(which(!x == ''))),
                  indicator = if(!is.null(cb$indicator)) as.character(cb$indicator[filter]) else '',
                  condition = if(!is.null(cb$condition)) as.character(cb$condition[filter]) else '')
  
  if(sum(duplicated(cb$code)) > 0) stop('Duplicate code labels not allowed:\n\t- ', paste(cb$code[duplicated(cb$code)], collapse='\n\t- '))
  
  cb$parent = NA
  parents = NA
  for(i in 1:nrow(cb)){
    j = cb$level[i]
    parents[j] = as.character(cb$code[i])
    cb$parent[i] = ifelse(j > 1, parents[j-1], NA)
  }
  
  cb = cb[,c('parent','code', 'indicator','condition')]
  class(cb) = c('codebook', class(cb))
  cb
}


plot.codebook <- function(cb, labelsize=0.5){
  rootvertices = cb$code[is.na(cb$parent)]
  vmeta = data.frame(name=unique(as.character(cb$code)))
  vmeta$root = ifelse(vmeta$name %in% rootvertices, T, F)
  g = graph.data.frame(cb[!is.na(cb$parent), c('parent','code')], vertices = vmeta)
  g$layout = layout_as_tree(g, root = which(V(g)$root), flip.y = F, circular=F)
  g$layout = cbind(g$layout[,2], g$layout[,1])

  plot(g, vertex.size=0, edge.arrow.size=0.00001, vertex.label.cex=labelsize)
}

codebookFromDict <- function(dict){
  if(!'dictionary' == class(dict)) stop('Dictionary has to have the quanteda dictionary class')
  codebook = sapply(dict, paste, collapse=' ')
  codebook = data.frame(parent=NA, code=names(codebook), indicator=as.character(codebook), condition='')
  class(codebook) = c('codebook', 'data.frame')
  codebook
}

