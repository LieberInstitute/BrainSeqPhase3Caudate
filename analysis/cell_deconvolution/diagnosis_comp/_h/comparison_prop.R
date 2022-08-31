## Calculate comparison between diagnosis per cell type
suppressPackageStartupMessages({
    library(tidyverse)
    library(ggpubr)
})

## Functions
save_img <- function(image, fn, w, h){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

get_pheno <- function(){
    fn <- paste0("/ceph/projects/v4_phase3_paper/inputs/",
                 "phenotypes/_m/merged_phenotypes.csv")
    return(data.table::fread(fn) %>%
           filter(Dx %in% c("SZ", "CTL"), Age > 17))
}
memPHENO <- memoise::memoise(get_pheno)

prep_data <- function(){
    load("../../_h/est_prop_Bisque.Rdata")
    return(est_prop_bisque$caudate$Est.prop.long %>%
           inner_join(memPHENO(), by=c("sample"="RNum")) %>%
           mutate_if(is.character, as.factor) %>%
           rename("Proportion"="prop"))
}
memDT <- memoise::memoise(prep_data)

cal_wilcox <- function(celltype){
    ctl <- memDT() %>%
        filter(cell_type == celltype, Dx == "CTL") %>%
        pull(Proportion)
    sz <- memDT() %>%
        filter(cell_type == celltype, Dx == "SZ") %>%
        pull(Proportion)
    res <- wilcox.test(ctl, sz, alternative="two.sided")
    return(list("Statistic"=res$statistic, "P_Value"=res$p.value))
}

#### MAIN
ct_lt = c(); stat_lt = c(); pval_lt = c();
for(celltype in unique(memDT()$cell_type)){
    l       <- cal_wilcox(celltype)
    ct_lt   <- c(ct_lt, celltype)
    stat_lt <- c(stat_lt, l[["Statistic"]])
    pval_lt <- c(pval_lt, l[["P_Value"]])
}
fn = "diagnosis_comparison_celltype_prop_2side_wilcox.txt"
data.frame("Celltype"=ct_lt, "Statistic"=stat_lt, "P_Value"=pval_lt) %>%
    data.table::fwrite(fn, sep='\t')

#### Reproducibility Information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
