\name{BIDCHIPS-package}
\alias{BIDCHIPS-package}
\alias{BIDCHIPS}
\docType{package}
\title{
BDCHIPS: Bias-Decomposition of ChIP-Seq, a reference implementation of the method of Perkins, Palidwor and Ramachandran for quantifying and removing biases from ChIP-seq data.
}
\description{
This package accompanies the paper "BIDCHIPS: Bias decomposition and removal from ChIP-seq data clarifies true binding signal and its functional correlates" by Ramachandran, Palidwor and Perkins, Epigenetics & Chromatin, Epigenetics & Chromatin 2015, 8:33. It uses a linear model to separate out non-binding influences from the ChIP-Seq signal based on such inputs as GC content, chromatin accessibility, IgG controls and input DNA.

}
\details{
\tabular{ll}{
Package: \tab BIDCHIPS\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2015-09-18\cr
License: \tab GPL-2\cr
}
}
\author{
Author: Gareth Palidwor
\cr
Maintainer: Gareth Palidwor <gpalidwor@ohri.ca>
}
\references{
"BIDCHIPS: Bias decomposition and removal from ChIP-seq data clarifies true binding signal and its functional correlates" by Ramachandran, Palidwor and Perkins, Epigenetics & Chromatin, Epigenetics & Chromatin 2015, 8:33
}
\keyword{ 
}
\seealso{
}
\examples{
#code for loading human hg19 chromosome 22 windows of data and generating linear model
#data downloadable from Perkins' Lab web site (http://dropbox.ogic.ca/Perkins_BIDCHIPS/hg19chr22BIDCHIPS_rawdata.zip)
#once the linear model and data are generated, we can apply the linear model to the data
#using the predict function.
chr22WindowModel<-
    buildModel(
      MappabilityBwFile="wgEncodeCrgMapabilityAlign36mer_chr22.bw",
      inBamFilesList=list(
        FNmDNaseI="wgEncodeUwDnaseGm12878AlnRep1_chr22.bam",
        FNmiDNA="wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22.bam",
        FNmIgGCtrl="wgEncodeSydhTfbsGm12878InputIggrabAlnRep1_chr22.bam"),
      outBamFileList=list(
        FNmTrt="wgEncodeHaibTfbsGm12878Atf3Pcr1xAlnRep1_chr22.bam"
      ),
      modeloutfile="hg19chr22model.object",
      windowsoutfile="hg19chr22windows.object",
      regions="chr22",ReadLen=36,cores=8,WinSz=257
)



#code for loading human hg19 windows of data based on peak regions and generating linear model and data windows
#data downloadable from Perkins' Lab web site (http://dropbox.ogic.ca/Perkins_BIDCHIPS/hg19allBIDCHIPS_rawdata.zip)
summits<-read.table("../data/wgEncodeHaibTfbsGm12878Atf3Pcr1xAlnRep1_summits.bed",header=F)[,1:3]
colnames(summits)<-c("PeakChr","SummitStart0based","SummitEnd1based")
summits<-subset(summits,!(PeakChr %in% c("chrM","chrY")))
summit_order<-order(match(summits$PeakChr,allchroms),summits$SummitStart0based)
summits<-summits[summit_order,]
peaksummit_gr<-GRanges(seqnames=summits$PeakChr,
                       IRanges(start=summits$SummitStart0based-HalfWin,end=summits$SummitStart0based+HalfWin),
                        seqlengths=chromlens[allchroms])
peakModel<-
  buildModel(
    MappabilityBwFile="wgEncodeCrgMapabilityAlign36mer.bw",
    inBamFilesList=list(
      FNmDNaseI="wgEncodeUwDnaseGm12878AlnRep1.bam",
      FNmiDNA="wgEncodeSydhTfbsGm12878InputStdAlnRep1.bam",
      FNmIgGCtrl="wgEncodeSydhTfbsGm12878InputIggrabAlnRep1.bam"),
    outBamFileList=list(
      FNmTrt="wgEncodeHaibTfbsGm12878Atf3Pcr1xAlnRep1.bam"
    ),
    modeloutfile="hg19peakmodel.object",
    windowsoutfile="hg19peakwindows.object",
    regions=peaksummit_gr,ReadLen=36,cores=2,WinSz=257
  )


#code for loading all human hg19 windows of data data and generating linear model and data windows
#data downloadable from Perkins' Lab web site (http://dropbox.ogic.ca/Perkins_BIDCHIPS/hg19chr22BIDCHIPS_rawdata.zip)
library(BSgenome.Hsapiens.UCSC.hg19)
allchroms<-grep("chr[123456789X0]+$",seqnames(BSgenome.Hsapiens.UCSC.hg19),value=T)
allWindowModel<-
  buildModel(
    MappabilityBwFile="wgEncodeCrgMapabilityAlign36mer.bw",
    inBamFilesList=list(
      FNmDNaseI="wgEncodeUwDnaseGm12878AlnRep1.bam",
      FNmiDNA="wgEncodeSydhTfbsGm12878InputStdAlnRep1.bam",
      FNmIgGCtrl="wgEncodeSydhTfbsGm12878InputIggrabAlnRep1.bam"),
    outBamFileList=list(
      FNmTrt="wgEncodeHaibTfbsGm12878Atf3Pcr1xAlnRep1.bam"
    ),
    modeloutfile="hg19allmodel.object",
    windowsoutfile="hg19allwindows.object",
    regions=allchroms,ReadLen=36,cores=2,WinSz=257
  )
predictedVals<-predict(chr22WindowModel$model,chr22WindowModel$windows)

#code for applying a model to windows and generating a bigwig file from it
#for visualization
library(stats)
predictedVals<-predict(chr22WindowModel$model,chr22WindowModel$windows)
residualVals<-chr22WindowModel$windows$FNmTrt-predictedVals
predicted_gr<-GRanges(seqnames=chr22WindowModel$windows$chrom,
                      IRanges(
                      start=chr22WindowModel$windows$start,
                      end=chr22WindowModel$windows$end),
                      score=predictedVals,
                      seqlengths=seqlengths(BSgenome.Hsapiens.UCSC.hg19))
residual_gr<-GRanges(seqnames=chr22WindowModel$windows$chrom,
                      IRanges(
                      start=chr22WindowModel$windows$start,
                      end=chr22WindowModel$windows$end),
                      score=residualVals,
                      seqlengths=seqlengths(BSgenome.Hsapiens.UCSC.hg19))
observed_gr<-GRanges(seqnames=chr22WindowModel$windows$chrom,
                      IRanges(
                      start=chr22WindowModel$windows$start,
                      end=chr22WindowModel$windows$end),
                      score=chr22WindowModel$windows$FNmTrt,
                      seqlengths=seqlengths(BSgenome.Hsapiens.UCSC.hg19))
#note that generate bigwig files can be quite big
export.bw(predicted_gr,"predicted.bw")
export.bw(observed_gr,"observed.bw")
export.bw(residual_gr,"residual.bw")
}