#' Generate a coverage object from a BAM file
#' 
#' @param bamfilename path to the BAM file
#' @return coverage object for whole BAM file
#' @examples
#' getBAMCoverage("/myproject/myBAM.bam")
getBAMCoverage<-function(bamfilename,regions){
  #add code so that it only returns full chromsomes for what we need
  what <- c("rname", "strand", "pos","qwidth")
  flag=scanBamFlag(isUnmappedQuery=FALSE,isSecondaryAlignem=FALSE)
  param<-ScanBamParam(what=what,flag)
  input_bam<-scanBam(bamfilename,param=param)
  b1<-lapply(names(input_bam[[1]]),function(elt){ do.call(c,unname(lapply(input_bam,"[[",elt)))})
  names(b1) <- names(input_bam[[1]])
  df2 <- do.call("DataFrame",b1)
  df2$rname<-factor(levels(input_bam[[1]]$rname)[df2$rname],levels=levels(input_bam[[1]]$rname)) 
  df2$strand<-factor(levels(input_bam[[1]]$strand)[df2$strand],levels=levels(input_bam[[1]]$strand))
  df2$qwidth<-input_bam[[1]]$qwidth
  #uniquified based on rname,strand,pos NOT qwidth
  #df2<-subset(df2,rname!="chrM")
  df2<-df2[!duplicated(df2[,1:3]),] 
  thestart=df2$pos+(df2$qwidth-1)*ifelse(df2$strand=="+",0,1)
  chromlens<-seqlengths(BSgenome.Hsapiens.UCSC.hg19)
  seqpos<-GRanges(seqnames=df2$rname,
                  IRanges(start=thestart,end=thestart),
                  seqlengths=chromlens[levels(df2$rname)])

  c1<-coverage(seqpos)
  c1
}  

#' Generate counts from a set of windows from a BAM file from a given set of chromosomes
#' 
#' @param bamfilename path to the BAM file
#' @param regions a list of chromosomes (seqnames) to provide window counts for
#' @param window size (window suze should be an odd integer value)
#' @param cores number of cores to parallelize using
#' @return list of read counts for windows in each chromosome 
#' getBAMCoverage("/myproject/myBAM.bam",c("chr22",257))
getBAMWindows<-function(bamfilename,regions,WinSz,cores){
  c1<-getBAMCoverage(bamfilename)
  rc<-list()
  for(chromname in regions){
    v1<-as.numeric(c1[[chromname]]) 
    rc[[chromname]]<-unlist(mclapply(1:floor(length(c1[[chromname]])/WinSz),
                                   function(x){sum(v1[((x-1)*WinSz+1):(x*WinSz)]/WinSz)},mc.cores=cores))
  }
  rc
}

#' Generate counts from a set of windows in a GRanges object
#' 
#' @param bamfilename path to the BAM file
#' @param GRanges object of windows
#' @return list of read counts for windows in each chromosome in order of the GRanges object
#' @examples
#' getBAMCoverage("/myproject/myBAM.bam",myGRanges)
getBAMGRanges<-function(bamfilename,peaksummits_gr,cores){
  c1<-getBAMCoverage(bamfilename)
  prc<-list()
  for(chromname in levels(seqnames(peaksummits_gr))){
    rr1<-as.data.frame(ranges(subset(peaksummits_gr,seqnames==chromname)))    
    v1<-as.numeric(c1[[chromname]]) 
    prc[[chromname]]<-unlist(mclapply(1:nrow(rr1),
      function(x){sum(v1[rr1[x,1]:rr1[x,2]])/(rr1[x,2]-rr1[x,1]+1)},mc.cores=cores))
  }
  prc
}

#' Return coverage object for mappability file
#' 
#' @param MappabilityBwFile path to the mappability BigWig File
#' @return mappability coverage object filted for mappability of "1" regions
#' @examples
#' getBAMCoverage("/myproject/myMappability.bw")
getMappabilityCoverage<-function(MappabilityBwFile){
  bwf<-BigWigFile(MappabilityBwFile)
  bw_rd<-import(bwf)  
  bw1_rd<-bw_rd[bw_rd$score==1]
  c1<-coverage(bw1_rd)
  c1
}


#' Return normalized mappability based on provided regions
#' 
#' @param MappabilityBwFile path to the mappability BigWig File
#' @param cores number of cores to parallelize using
#' @return mappability coverage object filtered for mappability of "1" regions
#' @examples
#' getBAMCoverage("/myproject/myMappability.bw",myGRanges,36)
getMappabilityGRanges<-function(MappabilityBwFile,regions,ReadLen,cores){
  c1<-getMappabilityCoverage(MappabilityBwFile)
  mappability<-list()
  allchromtable<-as.data.frame(regions)
  for(chromname in levels(allchromtable$seqnames)){
    print(chromname)
    chromtable<-subset(allchromtable,seqnames==chromname)
    chromcount<-mclapply(1:nrow(chromtable),
                         function(x){sum(as.numeric(c1[[chromname]][chromtable[x,2]:chromtable[x,3]]),
                                         as.numeric(c1[[chromname]][chromtable[x,2]:chromtable[x,3]+ReadLen-1])
                         )/(2*WinSz)},
                         mc.cores=cores)
    mappability[[chromname]]<-chromcount
  }
  mappability  
}

#' Return coverage object for mappability file based on provided regions
#' 
#' @param MappabilityBwFile path to the mappability BigWig File
#' @param cores number of cores to parallelize using
#' @return mappability coverage object filtered for mappability of "1" regions
#' @examples
#' getBAMCoverage("/myproject/myMappability.bw",myGRanges,36)
getMappabilityWindows<-function(MappabilityBwFile,regions,ReadLen,WinSz,cores){
  c1<-getMappabilityCoverage(MappabilityBwFile)
  mappability<-list()
  mpos_count<-list()
  mneg_count<-list()  
  for(chromname in regions){
    print(chromname)
    v1<-as.numeric(c1[[chromname]])
    mpos_count[[chromname]]<-mclapply(1:floor(length(c1[[chromname]])/WinSz),
                                      function(x){sum(v1[((x-1)*WinSz+1):(x*WinSz)])},
                                      mc.cores=cores)
    mneg_count[[chromname]]<-c(sum(v1[1:(WinSz-(ReadLen-1))]/(WinSz-(ReadLen-1))), #this is the first window
                               mclapply(2:floor(length(c1[[chromname]]-(ReadLen-1))/WinSz), 
                                        function(x){sum(v1[((x-1)*WinSz+1):(x*WinSz)-(ReadLen-1)])},
                                        mc.cores=cores))
    mappability[[chromname]]<-(unlist(mpos_count[[chromname]])+unlist(mneg_count[[chromname]]))/(2*WinSz)
  }
  mappability
}

#' Return list containing a model and the generated window values based on provided regions
#' 
#' @param MappabilityBwFile path to the mappability BigWig File
#' @param inBamFilesList a list of BAM files whose counts are to be used as the input to the linear model
#' @param outBamFilesList a list of one BAM file to be used as the output of the linear model
#' @param modeloutfile path to a file to contain the generated linear model
#' @param windowsoutfile path to a file to contain the generated windows
#' @param regions regions to analyze; either a list of chromosome names or a list of GenomicRanges 
#' @param ReadLen read length
#' @param window size should be an odd integer number. Only used if regions is a list of chromosomes
#' @return returns a list with the first element the generated model and the second element is the generated window values
#' @examples 
#' chr22WindowModel<-
#' buildModel(
#'   MappabilityBwFile="../data/wgEncodeCrgMappabilityAlign36mer_chr22.bw",
#'   inBamFilesList=list(
#'     FNmDNaseI="../data/wgEncodeUwDnaseGm12878AlnRep1_chr22.bam",
#'     FNmiDNA="../data/wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22.bam",
#'     FNmIgGCtrl="../data/wgEncodeSydhTfbsGm12878InputIggrabAlnRep1_chr22.bam"),
#'   outBamFileList=list(
#'     FNmTrt="../data/wgEncodeHaibTfbsGm12878Atf3Pcr1xAlnRep1_chr22.bam"
#'   ),
#'   modeloutfile="hg19chr22model.object",
#'   windowsoutfile="hg19chr22windows.object",
#'   regions="chr22",ReadLen=36,cores=2,WinSz=257
#' )
buildModel<-function(MappabilityBwFile,
                inBamFilesList,outBamFileList,
                modeloutfile,windowsoutfile,
                regions,ReadLen,cores,WinSz){
  
  require(BSgenome.Hsapiens.UCSC.hg19)
  require(GenomicRanges)
  require(rtracklayer)
  require(Rsamtools)
  require(biovizBase)
  require(parallel)
  
  if(class(regions) != "GRanges"){
    chromlens<-seqlengths(BSgenome.Hsapiens.UCSC.hg19)
    chroms<-regions
    allintstarts<-c()
    allintchrom<-c()
    for(chrom in chroms){
      chromlen<-chromlens[chrom]
      intstarts<-seq(1,chromlen,WinSz)
      allintstarts<-c(allintstarts,intstarts)
      intchrom<-rep(chrom,length(intstarts))
      allintchrom<-c(allintchrom,intchrom)
    }
    allintends<-allintstarts+WinSz-1
    g1<-GRanges(seqnames=allintchrom,
                IRanges(start=allintstarts,end=allintends),
                seqlengths=chromlens[chroms])
    g1<-trim(g1)
    g1<-g1[width(g1)==WinSz]
  }
  
  af1_gc<-c()
  positions<-data.frame()
  if(class(regions)=="GRanges"){
    af1_gc<-GCcontent(Hsapiens,regions)
    positions<-data.frame(chrom=seqnames(regions),start=start(regions),end=end(regions))
  }else{
    af1_gc<-GCcontent(Hsapiens, g1)
    positions<-data.frame(chrom=seqnames(g1),start=start(g1),end=end(g1))
  }

  inbamreadslist<-list()
  for(bam in names(inBamFilesList)){
    bamfile<-inBamFilesList[[bam]]
    if(class(regions)=="GRanges"){
      inbamreadslist[[bam]]<-getBAMGRanges(bamfile,regions,cores)
    }else{
      inbamreadslist[[bam]]<-getBAMWindows(bamfile,regions,WinSz,cores)
    }
  }
  
  outbamreadslist<-list()
  bam<-names(outBamFileList[1])
  bamfile<-outBamFileList[[bam]]
  if(class(regions)=="GRanges"){
    outbamreadslist[[bam]]<-getBAMGRanges(bamfile,regions,cores)
  }else{
    outbamreadslist[[bam]]<-getBAMWindows(bamfile,regions,WinSz,cores)
  }
  mappability<-c()
  if(class(regions)=="GRanges"){
    mappability<-getMappabilityGRanges(MappabilityBwFile,regions,ReadLen,cores)  
  }else{
    mappability<-getMappabilityWindows(MappabilityBwFile,regions,ReadLen,WinSz,cores)
  }
  df<-data.frame(mappability=unlist(mappability),gc=af1_gc[,1])
  df<-cbind(df,data.frame(lapply(inbamreadslist,function(x){unlist(x)})))
  df<-cbind(data.frame(lapply(outbamreadslist,function(x){unlist(x)})),df)
  formula<-paste(colnames(df)[1],"~",".")
  mymodel<-lm(formula,df)
  save(mymodel,file=modeloutfile)
  df<-data.frame(positions,df)
  row.names(df)<-NULL
  save(df,file=windowsoutfile)
  list(model=mymodel,windows=df)
}

