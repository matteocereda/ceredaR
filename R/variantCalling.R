# MUTECT 1.17 =====================

get_mutect = function( folder, pattern='.mutect.1.17.vcf'){
  f = list.files(folder, pattern, full.names = T)
  n = list.files(folder, pattern, full.names = F)
  n = gsub(pattern,"",n)
  names(f)=n
  hs =lapply(f, read.delim, comment.char = "#", h=T, sep="\t");names(hs)=n
  hs = lapply(hs, function(vvcf){
    vvcf$key=paste0(vvcf$contig,":",vvcf$position,"-",vvcf$position,"_",vvcf$ref_allele,"_",vvcf$alt_allele)
    return(vvcf)
  })
  hs
}

filter_mutect_1 = function(PWD, MIN_ALLELE_FREQ = 0.001, MIN_COV_POSITION = 10){
  print('Reading mutect results ...')
  mut = get_mutect(PWD)
  if(length(mut)>0){
    print('Filtering ...')
    k = lapply(mut, function(x) subset(x, judgement=="KEEP"
                                       | (judgement=="REJECT" & failure_reasons=='possible_contamination') ) )
    print('Strand bias check ...')
    s = lapply(k, function(x){
      sb = strsplit(x$strand_bias_counts,"[[:punct:]]")
      sb = lapply(sb, as.numeric)
      sb = lapply(sb, na.omit)
      sb = sapply(sb, function(x) (x[3]>0 & x[4]>0) )
      return(x[which(sb==T),])
    })
    print('VAF and COV checks ...')
    s = lapply(s, function(x){
      x$tot_cov = x$t_ref_count+x$t_alt_count
      return(subset(x,
                    tumor_f>=MIN_ALLELE_FREQ
                    &
                      tot_cov>MIN_COV_POSITION
      ))
    })
    return(s)
  }else{
    stop(message("Check your path"))
  }
}

# VARSCAN 2 =====================

read.varscan.vcf=function(vcffile, sample, exome.gr=NULL) {
  suppressPackageStartupMessages(require(VariantAnnotation))
  vvcf=readVcf(vcffile, genome="hg19")
  t.depth=geno(vvcf)$DP[,1]
  t.ref.count=geno(vvcf)$RD[,1]
  t.alt.count=geno(vvcf)$AD[,1]
  t.freq=as.numeric(gsub("%", "", geno(vvcf)$FREQ))/100

  ADF=geno(vvcf)$ADF[,1]
  ADR=geno(vvcf)$ADF[,1]

  v.info=as.data.frame(vvcf@rowRanges)
  REF=as.data.frame(vvcf@fixed$REF)[,1]
  ALT=as.data.frame(vvcf@fixed$ALT)[,3]
  ss = unlist(strsplit(sample,"\\."))[1]
  return(data.frame(v.info[,1:3], REF, ALT, DP= t.depth, DP.R=t.ref.count, DP.A=t.alt.count,t.freq
                    , ADF=ADF, ADR=ADR, sample=ss))
}

read.varscan.somatic.vcf=function(vcffile, sample, exome.gr=NULL) {
  suppressPackageStartupMessages(require(VariantAnnotation))
  ##fileformat=VCFv4.1
  ##source=VarScan2
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth of quality bases">
  ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description="Indicates if record is a somatic mutation">
  ##INFO=<ID=SS,Number=1,Type=String,Description="Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknown)">
  ##INFO=<ID=SSC,Number=1,Type=String,Description="Somatic score in Phred scale (0-255) derived from somatic p-value">
  ##INFO=<ID=GPV,Number=1,Type=Float,Description="Fisher's Exact Test P-value of tumor+normal versus no variant for Germline calls">
  ##INFO=<ID=SPV,Number=1,Type=Float,Description="Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls">
  ##FILTER=<ID=str10,Description="Less than 10% or more than 90% of variant supporting reads on one strand">
  ##FILTER=<ID=indelError,Description="Likely artifact due to indel reads at this position">
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
  ##FORMAT=<ID=RD,Number=1,Type=Integer,Description="Depth of reference-supporting bases (reads1)">
  ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
  ##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
  ##FORMAT=<ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev">
  vvcf=as.data.frame(readVcfAsVRanges(vcffile, genome="hg19"))
  colnames(vvcf)[1]='chr'
  vvcf$key=paste0(vvcf$chr,":",vvcf$start,"-",vvcf$end,"_",vvcf$ref,"_",vvcf$alt)

  vvcf$FREQ=as.numeric(gsub("%", "", vvcf$FREQ))/100

  t = subset(vvcf, sampleNames=="TUMOR")
  n = subset(vvcf, sampleNames=="NORMAL")

  t$totalDepth.NORMAL = n$totalDepth[match(t$key, n$key)]
  t$refDepth.NORMAL = n$refDepth[match(t$key, n$key)]
  t$altDepth.NORMAL= n$altDepth[match(t$key, n$key)]
  t=t[order(t$FREQ,decreasing = T),]
  t$sample=unlist(strsplit(sample,"\\."))[1]
  t
  }


get_VarScan = function( folder, pattern='.varscan2.vcf'){
  f = list.files(folder, pattern, full.names = T)
  n = list.files(folder,pattern, full.names = F)
  n = gsub(pattern,"",n)
  names(f)=n
  hs = mapply(read.varscan.vcf,f,n, SIMPLIFY = F );
  names(hs)=n
  hs
}

get_VarScan_somatic = function(folder, pattern='.varscan2.vcf'){
  f = list.files(folder, pattern, full.names = T)
  n = list.files(folder,pattern, full.names = F)
  n = gsub(pattern,"",n)
  names(f)=n
  hs = mapply(read.varscan.somatic.vcf,f,n, SIMPLIFY = F );
  names(hs)=n
  hs
}

filter_varscan_2 = function(PWD, PATTERN='varscan', MIN_ALLELE_FREQ = 0.001, MIN_COV_POSITION = 10 ){
  print("Reading VARSCAN2 results ...")
  v = get_VarScan_somatic(PWD, PATTERN)

  snp=grep('snp',names(v))
  indel=grep('indel',names(v))

  s = v[snp]
  i = v[indel]

  names(s) = sapply(s, function(x) x$sample[[1]])
  names(i) = sapply(i, function(x) x$sample[[1]])

  # cn = intersect(colnames(s[[1]]),colnames(i[[1]]))
  cn =c(
   "chr","start","end","strand","ref","alt","totalDepth","refDepth","altDepth"
   ,"sampleNames" ,"DP","SOMATIC"
   ,"SS","SSC","GPV", "SPV","GT","GQ"  , "RD","FREQ","DP4",
    "key","totalDepth.NORMAL","refDepth.NORMAL","altDepth.NORMAL","sample")

  s = lapply(s, function(x,y){ x[,y]}, y=cn)
  i = lapply(i, function(x,y){ x[,y]}, y=cn)

  varscan = mapply(rbind, s, i, SIMPLIFY = F)
  print("Keeping SOMATIC ...")
  varscan = lapply(varscan, function(x) subset(x, SOMATIC))

  print("Filter checks ...")

  varscan = lapply(varscan, function(x){
    subset(x,
           DP>=MIN_COV_POSITION
           # & (ADF>0 & ADR>0)
           & FREQ>=MIN_ALLELE_FREQ)
  })
  varscan
}


# MERGE MUTECT VARSCAN 2 =====================


merge = function(x,y,z=NULL){
  options(stringsAsFactors = F)
  colnames(x)=c('chr','start','end','ref','var','freq','cov','dbsnp_site')
  colnames(y)=c('chr','start','end','ref','var','freq','cov')
  if(!is.null(z))  colnames(z)=c('chr','start','end','ref','var','freq','cov')
  y$chr=as.character(y$chr)
  if(!is.null(z)) z$chr=as.character(z$chr)

  x$ID = with(x, paste0(chr, ".",end,".",ref,".", var))
  y$ID = with(y, paste0(chr, ".",end,".",ref,".", var))
  if(!is.null(z)) z$ID = with(z, paste0(chr, ".",end,".",ref,".", var))

  if(!is.null(z)){
    id = unique(c(x$ID,y$ID,z$ID))
  } else {
    id = unique(c(x$ID,y$ID))
  }
  df = data.frame( id=id, chr=NA, start=NA, end=NA, ref=NA, var=NA, freq=NA, cov=NA, MUTECT=NA, VARSCAN=NA )

  df[which(id%in%x$ID ), 2:8] = x[,1:7]
  df$dbsnp_site=NA
  df$dbsnp_site[which(id%in%x$ID )] = x$dbsnp_site

  tmp = which(df$id%in%y$ID & is.na(df$MUTECT) )
  tmp2 = df$id[tmp]
  df[ tmp, 2:8]  = y[match(tmp2, y$ID),1:7]

  if(!is.null(z)) df[which(id%in%z$ID ), 2:8] = z[,1:7]

  df$MUTECT = df$id %in% x$ID
  if(!is.null(z)) {
    df$VARSCAN = (df$id %in% y$ID) | (df$id %in% z$ID)
  } else {
    df$VARSCAN = (df$id %in% y$ID)
  }

  df
}

merge_mutect_varscan_somatic = function(m,v,show.plot=TRUE){

  m = m[,c("key","tot_cov","tumor_f")]
  colnames(m)=c("key","DP","FREQ")
  v = v[,c("key","DP","FREQ")]

  require(gplots)
  ven = venn(list('Mutect1.17'=m$key, "VarScan2"=v$key), show.plot = show.plot)
  df = rbind()
  if(!is.null(ven)){

    A = cbind()
    B = cbind()
    C = cbind()
    ids = attr(ven,"intersections")

    if(length(ids[['Mutect1.17']])>0) A=cbind( subset(m, key%in%ids[['Mutect1.17']]), method='Mutect1.17')
    if(length(ids[['Mutect1.17:VarScan2']])>0) B=cbind( subset(v, key%in%ids[['Mutect1.17:VarScan2']]), method='Mutect1.17:VarScan2')
    if(length(ids[['VarScan2']])>0) C=cbind( subset(v, key%in%ids[['VarScan2']]), method='VarScan2')

    if(!is.null(A)) if(nrow(A)>0) df=rbind(df, A)
    if(!is.null(B)) if(nrow(B)>0) df=rbind(df, B)
    if(!is.null(C)) if(nrow(C)>0) df=rbind(df, C)


  }
  df
}


