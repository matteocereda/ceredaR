pkgs = c('VariantAnnotation'
         ,'RColorBrewer'
         ,'plyr'
         ,'ggplot2'
         ,'GenomicRanges'
         ,"reshape"
         ,'reshape2'
         ,'grid'
         ,'dplyr'
         , 'ggrepel'
         ,'data.table')


lib = installed.packages()
installed=pkgs %in% rownames(lib)
not_installed = which(!installed)
if(length(not_installed)>0){
  for(i in not_installed) install.packages(pkgs[i])
}


for(i in pkgs) suppressPackageStartupMessages(library(i, character.only =T))

options(datatable.fread.datatable=FALSE
        , stringsAsFactors = F)

nonsilent=c('frameshift substitution','nonframeshift substitution','nonsynonymous SNV','splicing','stopgain','stoploss')
chroms = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
           "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
           "chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")




# common =====================================

mapply_pb = function(FUN, X, Y,  ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- mapply(wrapper, X, Y, ...)
  close(pb)
  res
}

lapply_pb = function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}


write.bed=function(...){
  write.table(..., row.names=F, col.names=F, quote=F, sep="\t")
}

