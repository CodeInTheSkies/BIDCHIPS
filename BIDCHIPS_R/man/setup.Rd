% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/BIDCHIPS.R
\name{setup}
\alias{setup}
\title{Return coverage object for mapability file based on provided regions}
\usage{
setup(MapabilityBwFile, inBamFilesList, outBamFileList, modeloutfile,
  windowsoutfile, regions, ReadLen, cores, WinSz)
}
\arguments{
\item{MapabilityBwFile}{path to the mapability BigWig File}

\item{inBamFilesList}{a list of BAM files whose counts are to be used as the input to the linear model}

\item{modeloutfile}{path to a file to contain the generated linear model}

\item{windowsoutfile}{path to a file to contain the generated windows}

\item{regions}{regions to analyze; if}

\item{ReadLen}{read length}

\item{outBamFilesList}{a list of one BAM file to be used as the output of the linear model}

\item{window}{size should be an odd integer number. Only used if regions is a list of chromosomes}
}
\value{
returns a list with the first element the generated model and the second element is the generated window counts
}
\description{
Return coverage object for mapability file based on provided regions
}
\examples{
chr22WindowModel<-
setup(
  MapabilityBwFile="../data/wgEncodeCrgMapabilityAlign36mer_chr22.bw",
  inBamFilesList=list(
    FNmDNaseI="../data/wgEncodeUwDnaseGm12878AlnRep1_chr22.bam",
    FNmiDNA="../data/wgEncodeSydhTfbsGm12878InputStdAlnRep1_chr22.bam",
    FNmIgGCtrl="../data/wgEncodeSydhTfbsGm12878InputIggrabAlnRep1_chr22.bam"),
  outBamFileList=list(
    FNmTrt="../data/wgEncodeHaibTfbsGm12878Atf3Pcr1xAlnRep1_chr22.bam"
  ),
  modeloutfile="hg19chr22model.object",
  windowsoutfile="hg19chr22windows.object",
  regions="chr22",ReadLen=36,cores=2,WinSz=257
)
}

