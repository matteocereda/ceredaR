prepare_ANNOVAR_infile=function(x, filename){
  k = lapply(x, function(x){
    tmp = strsplit(x$key, '[[:punct:]]')
    tmp = tmp[sapply(tmp, length)==5]
    do.call(rbind,tmp)
  })
  write.bed(unique(do.call(rbind,k)), file=filename)
}

set_damaging = function(x, med=10, high=20){
  x$damaging=NA
  x$damaging[which(as.numeric(x$CADD13_PHRED)>med)] = "Medium"
  x$damaging[which(as.numeric(x$CADD13_PHRED)>high)] = "High"

  # set the CADD score of indels to "-"
  id = which( x$ExonicFunc.refGene%in%c("frameshift substitution","nonframeshift substitution","splicing") )

  x$damaging[id] = "High"
  x

}

get_HGNC_isoform = function(x){
  data("hgnc")
  # , y="/Users/mcereda/Lavoro/HGNC_1612.Rdata"
  # if(file.exists(y)){
  #   load(y)
  # }else{
  #   stop(message("Missing HGNC file"))
  # }
  x$hgnc_canonical_refseq=NA
  x$alternative_refseq=NA

  suppressPackageStartupMessages(require(GenomicRanges))

  hgnc = hgnc[which(!is.na(hgnc$chromosome)),]

  a = with(x, GRanges( seqnames = Chr, IRanges( start = Start, end = End ) ) )
  b = with(hgnc,
           GRanges( seqnames = chromosome, IRanges( start = tstart, end = tend ) ))
  ov = findOverlaps(a,b)

  x$hgnc_refseq_accession = NA
  x$hgnc_refseq_accession[queryHits(ov)] = hgnc$refseq_accession[subjectHits(ov)]
  # x$hgnc_refseq_accession = sapply(strsplit(x$hgnc_refseq_accession, "\\|"),"[[",1)

  AA = strsplit(x$AAChange.refGene, "\\,")
  names(AA) = as.character(1:length(AA))
  nm = x$hgnc_refseq_accession
  names(nm) = as.character(1:length(nm))

  x$hgnc_canonical_refseq=
    unlist(
      mapply(function(x,y){
        if(!is.na(x[1]) & !is.na(y)){
          id =grep(y,x)
          x = ifelse(length(id) > 0, paste(x[id], collapse = "|"), NA)
          return(x)
        }else{
          return(NA)
        }
      }, AA, nm, SIMPLIFY = F)
    )
  id = which(is.na(x$hgnc_canonical_refseq))
  x$alternative_refseq[id]=x$AAChange.refGene[id]
  return(x)
}

set_genes_in_haloplex = function(x, fname="~/UV2_remote/epigenetics/Capture_kits/Agilent_Haloplex.Hematopoietic/04818-1466508813_Amplicons.bed"){
  suppressPackageStartupMessages(require(GenomicRanges))
  target = read.delim(file=fname, skip=2, h=F)
  a = with(x, GRanges( seqnames = Chr, IRanges( start = Start, end = End ) ) )
  b = GRanges( seqnames = target[,1], IRanges( start = target[,2], end = target[,3] ) )
  ov = findOverlaps(a,b)
  x$in_gene_panel=F
  x$in_gene_panel[ queryHits(ov) ] = T
  x
}

set_annotations_after_ANNOVAR <- function(annot, my_cancer_site, target_file_bed="~/UV2/epigenetics/Capture_kits/Agilent_Haloplex.Hematopoietic/04818-1466508813_Amplicons.bed") {
  data(cancerGenes)
  data("panel_drug")
  chroms = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8",
             "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
             "chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

  cancer_genes_leu = unique(subset(cancerGenes, cancer_site%in%my_cancer_site)$symbol)
  print("Setting cancer genes ...")
  annot$gene_type=geneInfo$cancer_type[match(annot$Gene.refGene, geneInfo$symbol)]
  annot$gene_type[which((annot$gene_type=='rst'))] = NA
  annot$cancer_gene_site=NA
  annot$cancer_gene_site[which(annot$Gene.refGene%in%cancer_genes_leu)] = my_cancer_site

  print("Setting functional annotations ...")
  annot$ExonicFunc.refGene[which(annot$Func.refGene=='splicing')] = 'splicing'
  annot$AAChange.refGene[which(annot$Func.refGene=='splicing')] = annot$GeneDetail.refGene[which(annot$Func.refGene=='splicing')]
  annot$AAChange.refGene[which(annot$Func.refGene=='UTR3')] = annot$GeneDetail.refGene[which(annot$Func.refGene=='UTR3')]
  annot$AAChange.refGene[which(annot$Func.refGene=='UTR5')] = annot$GeneDetail.refGene[which(annot$Func.refGene=='UTR5')]

  annot$Gene.refGene = sapply(strsplit(annot$Gene.refGene,"[[:punct:]]"), "[[", 1)

  annot = subset( annot, Chr%in%chroms)

  print("Setting damaging effect ...")
  annot = set_damaging(annot, med=5, high=20)
  print("Setting canonical annotation ...")
  annot = get_HGNC_isoform(annot)
  if(!is.null(target_file_bed)){
    print("Setting on target genes ...")
    annot = set_genes_in_haloplex(annot, target_file_bed)
  }
  annot$druggable=annot$Gene.refGene%in%panel_drug$gene_symbol
  annot
}


