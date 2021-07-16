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

assays(cov_gtex_caudate)$raw_counts <- raw_counts

## Save for later
save(cov_gtex_caudate, file = here("processed-data", "cov_gtex_caudate.Rdata"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

