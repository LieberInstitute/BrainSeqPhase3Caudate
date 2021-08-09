library("here")
library("recount3")
library("megadepth")
library("SummarizedExperiment")
library("rtracklayer")
library("sessioninfo")

## For the cache
dir.create(here("processed-data"), showWarnings = FALSE)
dir.create(here("processed-data/recount3-cache"), showWarnings = FALSE)
options("recount3_cache" = here("processed-data/recount3-cache"))

## Locate GTEx caudate BigWig files and metadata
################################################

## Locate human projects
human_projects <- available_projects()

## Subset to GTEx brain
proj_info <- subset(human_projects,
    project == "BRAIN" & file_source == "gtex")
proj_info
#      project organism file_source      project_home project_type n_samples
# 8688   BRAIN    human        gtex data_sources/gtex data_sources      2931

## Create an RSE object with a small annotation
rse_brain <- create_rse(proj_info, annotation = "ercc")

## Subset to caudate
table(rse_brain$gtex.smtsd)
#
#                          Brain - Amygdala
#                                       163
#  Brain - Anterior cingulate cortex (BA24)
#                                       201
#           Brain - Caudate (basal ganglia)
#                                       273
#             Brain - Cerebellar Hemisphere
#                                       250
#                        Brain - Cerebellum
#                                       285
#                            Brain - Cortex
#                                       286
#              Brain - Frontal Cortex (BA9)
#                                       224
#                       Brain - Hippocampus
#                                       220
#                      Brain - Hypothalamus
#                                       221
# Brain - Nucleus accumbens (basal ganglia)
#                                       262
#           Brain - Putamen (basal ganglia)
#                                       221
#        Brain - Spinal cord (cervical c-1)
#                                       171
#                  Brain - Substantia nigra
#                                       154
rse_caudate <-
    rse_brain[, rse_brain$gtex.smtsd == "Brain - Caudate (basal ganglia)"]


## Now load the expressed regions from
## https://github.com/LieberInstitute/qsva_brain/blob/master/process_bsp1_and_bsp3.R#L88
load(
    '/dcl01/lieber/RNAseq/Datasets/BrainSeq_hg38/Phase3/count_data/degradation_rse_phase3_caudate.rda',
    verbose = TRUE
)

## Start building our RSE object with the GTEx metadata and the ERs coordinates
cov_gtex_caudate <- SummarizedExperiment(colData = colData(rse_caudate),
    rowRanges = rowRanges(cov_rse_caudate))

## Extract the ranges into a BED file
export.bed(rowRanges(cov_rse_caudate),
    here("processed-data/qsva_ers.bed"))

dir.create(here("processed-data/megadepth-cache"), showWarnings = FALSE)
start <- Sys.time()
raw_counts <- sapply(seq_len(ncol(cov_gtex_caudate)), function(i) {
    ## Add a random sleep to avoid any limits from SciServer
    Sys.sleep(round(runif(1, min = 1, max = 5)))

    ## Extract the data with megadepth and keep the files around
    x <-
        get_coverage(
            cov_gtex_caudate$BigWigURL[i],
            op = "sum",
            annotation = here("processed-data/qsva_ers.bed"),
            prefix = here(
                "processed-data/megadepth-cache",
                paste0(colnames(cov_gtex_caudate)[i], ".bw.sum")
            )
        )
    x$score
})
end <- Sys.time()
rownames(raw_counts) <- rownames(cov_gtex_caudate)
colnames(raw_counts) <- colnames(cov_gtex_caudate)

## See how much it took
end - start
# Time difference of 58.20705 mins

assays(cov_gtex_caudate)$raw_counts <- raw_counts

## Inspect the final object
cov_gtex_caudate
# class: RangedSummarizedExperiment
# dim: 1000 273
# metadata(0):
# assays(1): raw_counts
# rownames(1000): chr12:109737294-109737352 chrX:73684564-73684674 ...
#   chr10:103648132-103648295 chr7:102649026-102650348
# rowData names(2): name score
# colnames(273): GTEX-1C6VS-0011-R5a-SM-7MXU9.1
#   GTEX-1CB4J-0011-R5a-SM-9WG5Y.1 ... GTEX-1IL2V-0011-R5a-SM-ARL78.1
#   GTEX-1IY9M-0011-R5a-SM-A9SLX.1
# colData names(198): rail_id external_id ... recount_seq_qc.errq
#   BigWigURL

## Counts are in very different scales
summary(as.vector(assays(cov_gtex_caudate)$raw_counts))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#    0        0      755    22866     8174 16627770
summary(as.vector(assays(cov_rse_caudate)$counts))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#  0.0       4.2      15.3     789.9      62.4 1658055.9

## Note that the initial qSVA code uses counts from
## https://github.com/LieberInstitute/qsva_brain/blob/master/process_bsp1_and_bsp3.R#L135
## which are normalized by SPEAQeasy (well, the older code) to
## 40 million 100 bp reads (per strand)
## https://github.com/LieberInstitute/RNAseq-pipeline/blob/master/sh/step5-coverage.sh#L120

## We can scale the recount3/megadepth counts to 40 million 100 bp reads
## also using the AUC as shown at
## http://research.libd.org/recount3/reference/transform_counts.html.
## Scale the counts using the AUC
assays(cov_gtex_caudate)$counts <- transform_counts(cov_gtex_caudate)

## Now the counts are in a more similar scale: there are some differences
## due to stranded bigwigs in SPEAQeasy vs unstranded in recount3 +
## well, BSP3 vs GTEx.
## But these are the counts you likely want to use later for qSVA using
## code similar to https://github.com/LieberInstitute/qsva_brain/blob/master/process_bsp1_and_bsp3.R#L128-L145
summary(as.vector(assays(cov_gtex_caudate)$counts))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0     0.0     5.0   137.2    51.0 85923.0


## Save for later
save(cov_gtex_caudate, file = here("processed-data", "cov_gtex_caudate.Rdata"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.1.0 Patched (2021-05-18 r80330)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2021-07-16
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
#  Biobase              * 2.52.0   2021-05-19 [2] Bioconductor
#  BiocFileCache          2.0.0    2021-05-19 [2] Bioconductor
#  BiocGenerics         * 0.38.0   2021-05-19 [2] Bioconductor
#  BiocIO                 1.2.0    2021-05-19 [2] Bioconductor
#  BiocParallel           1.26.1   2021-07-04 [2] Bioconductor
#  Biostrings             2.60.1   2021-06-06 [2] Bioconductor
#  bit                    4.0.4    2020-08-04 [2] CRAN (R 4.1.0)
#  bit64                  4.0.5    2020-08-30 [2] CRAN (R 4.1.0)
#  bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
#  blob                   1.2.1    2020-01-20 [2] CRAN (R 4.1.0)
#  cachem                 1.0.5    2021-05-15 [2] CRAN (R 4.1.0)
#  cli                    3.0.0    2021-06-30 [2] CRAN (R 4.1.0)
#  cmdfun                 1.0.2    2020-10-10 [1] CRAN (R 4.1.0)
#  colorout               1.2-2    2021-05-25 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
#  crayon                 1.4.1    2021-02-08 [2] CRAN (R 4.1.0)
#  curl                   4.3.2    2021-06-23 [2] CRAN (R 4.1.0)
#  data.table             1.14.0   2021-02-21 [2] CRAN (R 4.1.0)
#  DBI                    1.1.1    2021-01-15 [2] CRAN (R 4.1.0)
#  dbplyr                 2.1.1    2021-04-06 [2] CRAN (R 4.1.0)
#  DelayedArray           0.18.0   2021-05-19 [2] Bioconductor
#  desc                   1.3.0    2021-03-05 [2] CRAN (R 4.1.0)
#  digest                 0.6.27   2020-10-24 [2] CRAN (R 4.1.0)
#  dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
#  ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
#  fansi                  0.5.0    2021-05-25 [2] CRAN (R 4.1.0)
#  fastmap                1.1.0    2021-01-25 [2] CRAN (R 4.1.0)
#  filelock               1.0.2    2018-10-05 [2] CRAN (R 4.1.0)
#  fs                     1.5.0    2020-07-31 [2] CRAN (R 4.1.0)
#  generics               0.1.0    2020-10-31 [2] CRAN (R 4.1.0)
#  GenomeInfoDb         * 1.28.1   2021-07-01 [2] Bioconductor
#  GenomeInfoDbData       1.2.6    2021-05-11 [2] Bioconductor
#  GenomicAlignments      1.28.0   2021-05-19 [2] Bioconductor
#  GenomicRanges        * 1.44.0   2021-05-19 [2] Bioconductor
#  ggplot2                3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
#  glue                   1.4.2    2020-08-27 [2] CRAN (R 4.1.0)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
#  here                 * 1.0.1    2020-12-13 [1] CRAN (R 4.1.0)
#  hms                    1.1.0    2021-05-17 [2] CRAN (R 4.1.0)
#  htmltools              0.5.1.1  2021-01-22 [2] CRAN (R 4.1.0)
#  htmlwidgets            1.5.3    2020-12-10 [2] CRAN (R 4.1.0)
#  httpuv                 1.6.1    2021-05-07 [2] CRAN (R 4.1.0)
#  httr                   1.4.2    2020-07-20 [2] CRAN (R 4.1.0)
#  IRanges              * 2.26.0   2021-05-19 [2] Bioconductor
#  jsonlite               1.7.2    2020-12-09 [2] CRAN (R 4.1.0)
#  later                  1.2.0    2021-04-23 [2] CRAN (R 4.1.0)
#  lattice                0.20-44  2021-05-02 [3] CRAN (R 4.1.0)
#  lifecycle              1.0.0    2021-02-15 [2] CRAN (R 4.1.0)
#  magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
#  Matrix                 1.3-4    2021-06-01 [3] CRAN (R 4.1.0)
#  MatrixGenerics       * 1.4.0    2021-05-19 [2] Bioconductor
#  matrixStats          * 0.59.0   2021-06-01 [2] CRAN (R 4.1.0)
#  megadepth            * 1.2.0    2021-05-19 [1] Bioconductor
#  memoise                2.0.0    2021-01-26 [2] CRAN (R 4.1.0)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
#  pillar                 1.6.1    2021-05-16 [2] CRAN (R 4.1.0)
#  pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
#  pkgload                1.2.1    2021-04-06 [2] CRAN (R 4.1.0)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 4.1.0)
#  promises               1.2.0.1  2021-02-11 [2] CRAN (R 4.1.0)
#  purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
#  R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
#  R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
#  R.utils                2.10.1   2020-08-26 [2] CRAN (R 4.1.0)
#  R6                     2.5.0    2020-10-28 [2] CRAN (R 4.1.0)
#  rappdirs               0.3.3    2021-01-31 [2] CRAN (R 4.1.0)
#  Rcpp                   1.0.7    2021-07-07 [2] CRAN (R 4.1.0)
#  RCurl                  1.98-1.3 2021-03-16 [2] CRAN (R 4.1.0)
#  readr                  1.4.0    2020-10-05 [2] CRAN (R 4.1.0)
#  recount3             * 1.2.1    2021-05-25 [1] Bioconductor
#  restfulr               0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
#  rjson                  0.2.20   2018-06-08 [2] CRAN (R 4.1.0)
#  rlang                  0.4.11   2021-04-30 [2] CRAN (R 4.1.0)
#  rmote                  0.3.4    2021-05-25 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
#  Rsamtools              2.8.0    2021-05-19 [2] Bioconductor
#  RSQLite                2.2.7    2021-04-22 [2] CRAN (R 4.1.0)
#  rstudioapi             0.13     2020-11-12 [2] CRAN (R 4.1.0)
#  rtracklayer          * 1.52.0   2021-05-19 [2] Bioconductor
#  S4Vectors            * 0.30.0   2021-05-19 [2] Bioconductor
#  scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
#  servr                  0.22     2021-04-14 [1] CRAN (R 4.1.0)
#  sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 4.1.0)
#  SummarizedExperiment * 1.22.0   2021-05-19 [2] Bioconductor
#  testthat               3.0.4    2021-07-01 [2] CRAN (R 4.1.0)
#  tibble                 3.1.2    2021-05-16 [2] CRAN (R 4.1.0)
#  tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
#  utf8                   1.2.1    2021-03-12 [2] CRAN (R 4.1.0)
#  vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
#  withr                  2.4.2    2021-04-18 [2] CRAN (R 4.1.0)
#  xfun                   0.24     2021-06-15 [2] CRAN (R 4.1.0)
#  XML                    3.99-0.6 2021-03-16 [2] CRAN (R 4.1.0)
#  XVector                0.32.0   2021-05-19 [2] Bioconductor
#  yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
#  zlibbioc               1.38.0   2021-05-19 [2] Bioconductor
#
# [1] /users/lcollado/R/4.1
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1/R/4.1/lib64/R/library
