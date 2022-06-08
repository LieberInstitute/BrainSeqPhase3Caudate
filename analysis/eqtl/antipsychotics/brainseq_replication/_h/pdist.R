## This script examines the p-value distributions of random
## DEGs that are associated only with SZ AP dysregulation compared
## with SZ no AP dysregulation. Replication in BrainSEQ Phase 2.
library(tidyverse)

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_eqtl <- function(fn){
    eqtl <- data.table::fread(fn)
    return(eqtl)
}
memEQTL <- memoise::memoise(load_eqtl)

load_topVariant <- function(fn){
    eqtl <- memEQTL(fn) %>% group_by(gene_id) %>% arrange(pval_nominal) %>%
        filter(row_number()==1) %>% as.data.frame
    return(eqtl)
}

load_degs <- function(label){
    lt = list("AP"=paste0("../../../../differential_expression/antipsychotics/",
                          "_m/genes/diffExpr_sz_APVctl_full.txt"),
              "noAP"=paste0("../../../../differential_expression/antipsychotics/",
                            "_m/genes/diffExpr_sz_noAPVctl_full.txt"),
              "DEG"=paste0("../../../../differential_expression/_m/genes/",
                           "diffExpr_szVctl_full.txt"))
    df <- data.table::fread(lt[[label]]) %>% filter(adj.P.Val < 0.05)
    return(df)
}

get_unique_sets <- function(){
                                        # Get DEGs
    df1 <- load_degs("AP")$gencodeID
    df2 <- load_degs("noAP")$gencodeID
    df0 <- load_degs("DEG")$gencodeID
                                        # Select unique to AP or no AP
    ap <- setdiff(setdiff(df1, df2), df0)
    noap <- setdiff(setdiff(df2, df1), df0)
    return(list("AP"=ap, "noAP"=noap))
}

dist_comp <- function(fn, label){
    ### General t test comparing p-values
    eqtl_df <- memEQTL(fn)
    deg_lt <- get_unique_sets()
    df1 <- eqtl_df %>% filter(gene_id %in% deg_lt[["AP"]]) %>%
        mutate(Status="AP", log10_pval=-log10(pval_nominal))
    df2 <- eqtl_df %>% filter(gene_id %in% deg_lt[["noAP"]]) %>%
        mutate(Status="No AP", log10_pval=-log10(pval_nominal))
                                        # Statistical testing
    ## Test to see if AP has larger p-values than no AP
    print(wilcox.test(df1$pval_nominal, df2$pval_nominal,
                      alternative="greater"))
                                        # Plot density of all genes
    df <- bind_rows(df1, df2) %>% select(Status, pval_nominal, log10_pval) %>%
        mutate_if(is.character, as.factor)
    dist <- ggpubr::ggdensity(df, x="log10_pval", fill="Status", add="mean",
                              rug=TRUE, palette="npg", color="Status",
                              xlab="-log10(P)",
                              ggtheme=ggpubr::theme_pubr(base_size=15))
    outfile <- paste0("density_plot_",label)
    save_ggplots(outfile, dist, 7, 4)

    ### Randomly compare genes
    set.seed(13131313)
    df1 <- eqtl_df %>% filter(gene_id %in% sample(deg_lt[["AP"]], 10)) %>%
        mutate(Status="AP", log10_pval=-log10(pval_nominal))
    df2 <- eqtl_df %>% filter(gene_id %in% sample(deg_lt[["noAP"]], 10)) %>%
        mutate(Status="No AP", log10_pval=-log10(pval_nominal))
                                        # Statistical testing
    ## Test to see if AP has larger p-values than no AP
    print(wilcox.test(df1$pval_nominal, df2$pval_nominal,
                      alternative="greater"))
    df <- bind_rows(df1, df2) %>% select(Status, pval_nominal, log10_pval) %>%
        mutate_if(is.character, as.factor)
                                        # Plot density of randomly selected genes
    dist <- ggpubr::ggdensity(df, x="log10_pval", fill="Status", add="mean",
                              rug=TRUE, palette="npg", color="Status",
                              xlab="-log10(P)",
                              ggtheme=ggpubr::theme_pubr(base_size=15))
    outfile <- paste0("density_plot_random_",label)
    save_ggplots(outfile, dist, 7, 4)
}

dist_comp_topVariant <- function(fn, label){
    ### General t test comparing p-values
    eqtl_df <- load_topVariant(fn)
    deg_lt <- get_unique_sets()
    df1 <- eqtl_df %>% filter(gene_id %in% deg_lt[["AP"]]) %>%
        mutate(Status="AP", log10_pval=-log10(pval_nominal))
    df2 <- eqtl_df %>% filter(gene_id %in% deg_lt[["noAP"]]) %>%
        mutate(Status="No AP", log10_pval=-log10(pval_nominal))
                                        # Statistical testing
    ## Test to see if AP has larger p-values than no AP
    print(wilcox.test(df1$pval_nominal, df2$pval_nominal,
                      alternative="greater"))
                                        # Plot density of all genes
    df <- bind_rows(df1, df2) %>% select(Status, pval_nominal, log10_pval) %>%
        mutate_if(is.character, as.factor)
    dist <- ggpubr::ggdensity(df, x="log10_pval", fill="Status", add="mean",
                              rug=TRUE, palette="npg", color="Status",
                              xlab="-log10(P)",
                              ggtheme=ggpubr::theme_pubr(base_size=15))
    outfile <- paste0("density_plot_topVariant_",label)
    save_ggplots(outfile, dist, 7, 4)
}

## BrainSEQ DLPFC
fn1 = paste0("/ceph/projects/v4_phase3_paper/analysis/eqtl_analysis/bsp2/",
             "dlpfc/genes/prepare_expression/fastqtl_nominal/_m/",
             "Brainseq_LIBD.allpairs.txt.gz")
label1 = "BS_DLPFC"
dist_comp(fn1, label1)
dist_comp_topVariant(fn1, label1)
## BrainSEQ Hippocampus
fn2 = paste0("/ceph/projects/v4_phase3_paper/analysis/eqtl_analysis/bsp2/",
             "hippocampus/genes/prepare_expression/fastqtl_nominal/_m/",
             "Brainseq_LIBD.allpairs.txt.gz")
label2 = "BS_Hippocampus"
dist_comp(fn2, label2)
dist_comp_topVariant(fn2, label2)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
