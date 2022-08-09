This folder contains R and Matlab codes in support of the paper:

[1] P.Ramachandran,G.Palidwor,T.J.Perkins(2015)“BIDCHIPS:biasdecompositionandremovalfrom ChIP-seq data clarifies true binding signal and its functional correlates.” Epigenetics & Chromatin, Vol 8, Art. 33 (16 pages).

If you use this code, please cite the paper above!

The Matlab code comes as a single m-file that has all the required dependencies within it. Detailed header comments are given at the beginning of the file explaining usage with examples.

The R code can be found in the folder BIDCHIPS_R.

All input ChIP-seq and other short-read sequencing data has to be in the BED format (https://genome.ucsc.edu/FAQ/FAQformat.html#format1). BED can easily be obtained from BAM files using the bam2bed utility: http://bedtools.readthedocs.org/en/latest/content/tools/bamtobed.html

The mappability and GC-content tracks can be obtained using the mappability intervals files (downloadable from the UCSC Table Browser), and the genomic sequence (also downloadable from the UCSC website).

The folder also includes a test datasets for human genome chromosome 22 (which is smaller than a whole genome dataset and runs quicker). Its files show you want input data should look like, but you don't need them to run the code if you have your own data.