## This script generates pretty venn diagrams with ggvenn.
suppressPackageStartupMessages({
    library(dplyr)
    library(ggvenn)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

get_mouse2human <- function(){
    fn = paste0("/ceph/users/jbenja13/projects/aanri/racial_diff/input/",
                "celltypes/_h/cell_type/mouse2human_gene.txt")
    return(data.table::fread(fn))
}
memMOUSE <- memoise::memoise(get_mouse2human)

get_brainseq <- function(){
    ap <- data.table::fread("../../../_m/genes/diffExpr_sz_APVctl_full.txt") %>%
        filter(adj.P.Val < 0.05)
    noap <- data.table::fread("../../../_m/genes/diffExpr_sz_noAPVctl_full.txt") %>%
        filter(adj.P.Val < 0.05)
    degs <- data.table::fread("../../../../_m/genes/diffExpr_szVctl_full.txt") %>%
        filter(adj.P.Val < 0.05)
    return(list("AP"=ap, "noAP"=noap, "SZ"=degs))
}

animal_studies <- function(){
    chong <- data.table::fread("../../_m/Chong2002_rat.csv") %>%
        rename("Symbol_human"="Gene Name (Human)")
    korostynski <- data.table::fread("../../_m/Korostynski2013_mice.csv") %>%
        inner_join(memMOUSE(), by=c("Gene_symbol"="Symbol_mouse")) %>%
        filter(ANOVA_drug_FDR < 0.05)
    kim <- data.table::fread("../../_m/Kim2018_mice.csv") %>%
        inner_join(memMOUSE(), by=c("Gene symbol"="Symbol_mouse")) %>%
        filter(Tissue == "Striatum")
    return(list("Chong 2002"=chong, "Korostynski 2013"=korostynski,
                "Kim 2018"=kim))
}

venn_diagrams <- function(lab1, lab2){
    outfile = paste("antipsychotics_venn", gsub(" ", "_", lab1),
                    gsub(" ", "_", lab2), sep="_")
    x = list(
        A = animal_studies()[[lab1]] %>% select(Symbol_human) %>% unlist(),
        B = get_brainseq()[[lab2]] %>% select(Symbol) %>% unlist()
    )
    names(x) <- c(lab1, paste0("BrainSEQ (",lab2,")"))
    vv <- ggvenn(x, fill_color=ggpubr::get_palette(palette="npg", 2),
                 stroke_size = 0.5)
    save_ggplots(tolower(outfile), vv, 5, 5)
}

## Generate pretty venn diagrams
for(lab1 in c("Chong 2002", "Korostynski 2013", "Kim 2018")){
    for(lab2 in c("SZ", "AP", "noAP")){
        venn_diagrams(lab1, lab2)
    }
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
