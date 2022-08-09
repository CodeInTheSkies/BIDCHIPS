## BIDCHIPS
 
#### BIDCHIPS: bias decomposition and removal from ChIP-seq data

### Citation
 
This folder contains R and Matlab codes in support of the publication:

[BIDCHIPS: bias decomposition and removal from ChIP-seq data clarifies true binding signal and its functional correlates](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-015-0028-2), _Epigenetics & Chromatin_, 2015.

If you use this code, please cite the above publication.

### Matlab version

[The Matlab code comes as a single m-file](/bidchips.m) that has all the required dependencies within it. Detailed header comments are given at the beginning of the file explaining usage with examples.

### R version

The R code can be found in the folder [BIDCHIPS_R](/BIDCHIPS_R), and was written and prepared by the work's co-author Gareth Palidwor.

All input ChIP-seq and other short-read sequencing data has to be in the [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). BED can easily be obtained from BAM files using the [bam2bed utility](http://bedtools.readthedocs.org/en/latest/content/tools/bamtobed.html).

The mappability and GC-content tracks can be obtained using the mappability intervals files (downloadable from the [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables)), and the genomic sequence (also downloadable from [UCSC](https://genome.ucsc.edu/FAQ/FAQdownloads.html)).

The folder also includes a [test dataset for human genome chromosome 22](/hg19chr22BIDCHIPS_rawdata) (which is smaller than a whole genome dataset and runs quicker). Its files show you what input data should look like, but you don't need them to run the code if you have your own data.
