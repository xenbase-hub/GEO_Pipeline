# RNA-SEQ ANALYSIS

## Prerequisites

- samtools/1.8.0
- deeptools/3.1.3
- bowtie2/2.1.0
- rsem/1.3.0
- R/4.1.1 (or later)
- R packages RUVSeq 1.29.0 and tximport 1.22.0

All reference files:
https://drive.google.com/file/d/1gtfdp3bDVWutqaeUxj-GbMhUp7hjiX2K/view?usp=sharing

# Usage

All FASTQ downloading/processing, alignment, bam processing, and peak calling should performed on the command line.

Bash scripts should be submitted through the scheduler of your choice.

## R session info
```
 package              * version  date (UTC) lib source
 AnnotationDbi          1.56.2   2021-11-09 [1] Bioconductor
 aroma.light            3.24.0   2021-10-26 [1] Bioconductor
 assertthat             0.2.1    2019-03-21 [1] CRAN (R 4.1.1)
 Biobase              * 2.54.0   2021-10-26 [1] Bioconductor
 BiocFileCache          2.2.0    2021-10-26 [1] Bioconductor
 BiocGenerics         * 0.40.0   2021-10-26 [1] Bioconductor
 BiocIO                 1.4.0    2021-10-26 [1] Bioconductor
 BiocParallel         * 1.28.3   2021-12-09 [1] Bioconductor
 biomaRt                2.50.3   2022-02-03 [1] Bioconductor
 Biostrings           * 2.62.0   2021-10-26 [1] Bioconductor
 bit                    4.0.4    2020-08-04 [1] CRAN (R 4.1.1)
 bit64                  4.0.5    2020-08-30 [1] CRAN (R 4.1.1)
 bitops                 1.0-7    2021-04-24 [1] CRAN (R 4.1.1)
 blob                   1.2.2    2021-07-23 [1] CRAN (R 4.1.1)
 cachem                 1.0.6    2021-08-19 [1] CRAN (R 4.1.1)
 callr                  3.7.0    2021-04-20 [1] CRAN (R 4.1.1)
 cli                    3.6.1    2023-03-23 [1] CRAN (R 4.1.1)
 crayon                 1.4.2    2021-10-29 [1] CRAN (R 4.1.1)
 curl                 * 4.3.2    2021-06-23 [1] CRAN (R 4.1.1)
 DBI                    1.1.1    2021-01-15 [1] CRAN (R 4.1.1)
 dbplyr                 2.1.1    2021-04-06 [1] CRAN (R 4.1.1)
 DelayedArray           0.20.0   2021-10-26 [1] Bioconductor
 desc                   1.4.0    2021-09-28 [1] CRAN (R 4.1.1)
 devtools             * 2.4.3    2021-11-30 [1] CRAN (R 4.1.1)
 digest                 0.6.31   2022-12-11 [1] CRAN (R 4.1.1)
 dplyr                  1.0.8    2022-02-08 [1] CRAN (R 4.1.1)
 EDASeq               * 2.28.0   2021-10-26 [1] Bioconductor
 edgeR                * 3.36.0   2021-10-26 [1] Bioconductor
 ellipsis               0.3.2    2021-04-29 [1] CRAN (R 4.1.1)
 fansi                  0.5.0    2021-05-25 [1] CRAN (R 4.1.1)
 fastmap                1.1.0    2021-01-25 [1] CRAN (R 4.1.1)
 filelock               1.0.2    2018-10-05 [1] CRAN (R 4.1.1)
 fs                     1.5.2    2021-12-08 [1] CRAN (R 4.1.1)
 generics               0.1.3    2022-07-05 [1] CRAN (R 4.1.1)
 GenomeInfoDb         * 1.30.0   2021-10-26 [1] Bioconductor
 GenomeInfoDbData       1.2.7    2021-12-20 [1] Bioconductor
 GenomicAlignments    * 1.30.0   2021-10-26 [1] Bioconductor
 GenomicFeatures        1.46.1   2021-10-27 [1] Bioconductor
 GenomicRanges        * 1.46.1   2021-11-18 [1] Bioconductor


