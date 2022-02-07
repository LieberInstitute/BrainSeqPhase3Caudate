## This script examines the p-value distributions of random
## DEGs that are associated only with SZ AP dysregulation compared
## with SZ no AP dysregulation. The p-value distributions are
## from the SNPxDiagnosis interaction FastQTL analysis.
library(tidyverse)

load_eqtl <- function(fn){
    eqtl <- data.table::fread(fn)
    return(eqtl)
}
memEQTL <- memoise::memoise(load_eqtl)

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

extract_eQTL <- function(fn){
    ### Load and extract eQTL (q-value < 0.05)
    eqtl_df <- memEQTL(fn) %>% filter(qval < 0.05)
    deg_lt <- get_unique_sets()
    df1 <- eqtl_df %>% filter(gene_id %in% deg_lt[["AP"]]) %>%
        mutate(Status="AP")
    df2 <- eqtl_df %>% filter(gene_id %in% deg_lt[["noAP"]]) %>%
        mutate(Status="No AP")
    df <- bind_rows(df1, df2) %>% mutate_if(is.character, as.factor)
    data.table::fwrite(df, "antipsychotics_annotated_eQTL.tsv", sep='\t')
}

#### Main section
fn = paste0("/ceph/projects/v4_phase3_paper/analysis/eqtl_analysis/",
            "control_only/genes/expression_gct/prepare_expression/",
            "fastqtl_permutation/_m/Brainseq_LIBD.genes.txt.gz")
extract_eQTL(fn)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
