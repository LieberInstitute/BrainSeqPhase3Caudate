suppressPackageStartupMessages({
    library(argparse)
    library(doParallel)
    library(foreach)
    library(mashr)
    library(dplyr)
})


get_pvals <- function(feature){
    fn = paste0("../../_m/", feature, "/pvalue_fastqtl_ancestry.tsv")
    cClasses = c("character", "character", "numeric", "numeric", "numeric")
    pval <- data.table::fread(fn, header=TRUE, sep="\t",
                              drop="ALL", colClasses=cClasses) %>%
        mutate(effect=paste(gene_id, variant_id, sep="_")) %>%
        group_by(gene_id) %>%
        arrange(AA, EA, .by_group = TRUE) %>% slice(1)
    return(pval)
}

get_bhat <- function(feature){
    fn = paste0("../../_m/", feature, "/bhat_fastqtl_ancestry.tsv")
    cClasses = c("character", "character", "numeric", "numeric", "numeric")
    return(data.table::fread(fn, header=TRUE, sep="\t",
                             drop="ALL", colClasses=cClasses) %>%
           mutate(effect=paste(gene_id, variant_id, sep="_")) %>%
           distinct(effect, .keep_all=TRUE))
}

get_shat <- function(feature){
    fn = paste0("../../_m/", feature, "/shat_fastqtl_ancestry.tsv")
    cClasses = c("character", "character", "numeric", "numeric", "numeric")
    return(data.table::fread(fn, header=TRUE, sep="\t",
                             drop="ALL", colClasses=cClasses) %>%
           mutate(effect=paste(gene_id, variant_id, sep="_")) %>%
           distinct(effect, .keep_all=TRUE))
}

get_rand_set <- function(bhat0, percentage, SEED){
    set.seed(SEED)
    rand_n = round(dim(bhat0)[1] * percentage)
    return(sample(bhat0$effect, rand_n))
}

save_results_chunk <- function(feature, chunk_size, threads){
    SEED = 13131313; percentage = 0.01
    ## Load prepared data
    bhat0 <- get_bhat(feature)
    shat0 <- get_shat(feature)
    pval <- get_pvals(feature)
    ## Get random and strong set
    strong_set <- pval$effect
    random_set <- get_rand_set(bhat0, percentage, SEED)
    rm(pval)
    ## Prepared data for mashr
    bhat <- bhat0 %>% tibble::column_to_rownames("effect") %>%
        select(-gene_id, -variant_id) %>% as.matrix
    shat <- shat0 %>% tibble::column_to_rownames("effect") %>%
        select(-gene_id, -variant_id) %>% as.matrix
    rm(bhat0, shat0)
    ## Calculate correlation matrix
    data.temp = mash_set_data(bhat[random_set,],shat[random_set,])
    Vhat = estimate_null_correlation_simple(data.temp)
    rm(data.temp)
    ## Initial mash
    data.random = mash_set_data(bhat[random_set,],shat[random_set,], V=Vhat)
    data.strong = mash_set_data(bhat[strong_set,],shat[strong_set,], V=Vhat)
    ## Calculate data driven covariances
    U.pca = cov_pca(data.strong, 2) # based on N conditions
    U.ed = cov_ed(data.strong, U.pca)
    U.c = cov_canonical(data.random)
    ## Fit mash model
    m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
    rm(U.pca, U.ed, U.c)
    ## Chunking
    numCores <- threads #detectCores()
    registerDoParallel(numCores)
    fn = paste0(feature,"/lfsr_allpairs_ancestry.txt.gz")
    chunks = split(1:dim(bhat)[1],
                   cut(1:dim(bhat)[1], chunk_size, labels=FALSE))
    dt <- foreach(i=1:chunk_size, .combine=rbind) %dopar% {
        data.chunk = mash_set_data(bhat[chunks[[i]],],
                                   shat[chunks[[i]],], V=Vhat)
        m.chunk = mash(data.chunk, g=get_fitted_g(m), fixg=TRUE)
        data.frame(m.chunk$result$lfsr) %>%
            tibble::rownames_to_column("effect") %>%
            tidyr::separate(effect, c("gene_id", "variant_id"),
                            remove=FALSE, sep="_")
    }
    dt %>% data.table::fwrite(fn, sep='\t')
}

## Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-r", "--run_chunk", action="store_true", default=FALSE,
                    help="Run this as a chunk [default]")
parser$add_argument("-c", "--chunk_size", type="integer", default=250,
                    help="Chunk size used for parallel run [default: %default]")
parser$add_argument("-t", "--threads", type="integer", default=10,
                    help="Number of threads to run on [default: %default]")
args <- parser$parse_args()

## Run mashr for specific feature
if(args$run_chunk){
    save_results_chunk(args$feature, args$chunk_size, args$threads)
    ## print(paste("Run chunk:", args$chunk))
}

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
