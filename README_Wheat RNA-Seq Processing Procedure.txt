The raw RNA-seq reads (available at GEO GSE137895) were preprocessed by trimming the adaptor sequences, filtering low-quality reads (Phred Score ≤ 20) and eliminating short reads (length ≤ 20 bps) using software package FASTX-toolkit.
The cleaned RNA-seq reads in each sample were mapped using STAR to generate gene-level counts.
Then run the R code (provided in github) for differential expression analysis using DEseq2.
