1. The raw RNA-seq reads (available at GEO GSE137895) were preprocessed by trimming the adaptor sequences, filtering low-quality reads (Phred Score ≤ 20 [78]) and eliminating short reads (length ≤ 20 bps) using software package FASTX-toolkit.
2. The cleaned RNA-seq reads in each sample were mapped using STAR to generate gene-level counts.
3. Then run the R code (provided in github) for differential expression analysis using DEseq2.

