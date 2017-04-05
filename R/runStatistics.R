# READ GATK files ============================
setdec <- function(lines) {
  res = "."
  if (length(grep(",", lines)) > 0) {
    res=","
  }
  res
}

hs_metrics = function(infile){
  lines <- readLines(infile)
  dec = setdec(lines[7:8])
  metrics <- rbind(read.table(textConnection(lines[7:8]), header=TRUE, fill=TRUE, sep="\t", na.strings=c("?"), dec=dec))
  as.data.frame(metrics, stringsAsFactors=F)
}

get_hs_metrics <- function(folder, sname="E") {
  f = list.files(folder, "hs_metrics", full.names = T)
  n = list.files(folder, "hs_metrics", full.names = F)
  n = gsub(".hs_metrics.txt","",n)
  names(f)=n
  f =f[grep(sname, f)]
  hs =lapply(f, hs_metrics);
  hs = do.call(rbind, hs)
  # hs$SAMPLE=rownames(hs)
  # hs$LIBRARY=date
  hs
}


get_targeted_metrics <- function(folder, sname="E") {
  f = list.files(folder, "targeted_metrics", full.names = T)
  n = list.files(folder, "targeted_metrics", full.names = F)
  n = gsub(".targeted_metrics.txt","",n)
  names(f)=n
  f =f[grep(sname, f)]
  hs =lapply(f, hs_metrics);

  hs = do.call(rbind, hs)
  # hs$SAMPLE=rownames(hs)
  # hs$LIBRARY=date
  hs
}


get_coverage = function(folder, sname="E"){
  f = list.files(folder, "coverage", full.names = T)
  n = list.files(folder, "coverage", full.names = F)
  n = gsub(".coverage.txt","",n)
  names(f) = n
  f =f[grep(sname, f)]
  hs =lapply(f, read.delim)#;names(hs) =n
  hs
}

