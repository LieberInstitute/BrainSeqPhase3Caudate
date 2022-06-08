## This script examines tissue specific genes for correlation
## with gene expression
library(magrittr)

get_eGenes <- function(){
    fn = "../../summary_table/_m/brain_region_specific_eGenes.tsv"
    return(data.table::fread(fn))
}

get_pheno <- function(){
    fn = "/ceph/projects/v4_phase3_paper/inputs/phenotypes/_m/merged_phenotypes.csv"
    df = data.table::fread(fn) %>%
        dplyr::filter(Age > 13, Race %in% c("AA", "EA"), Dx %in% c("CTL", "SZ"))
    return(df)
}
memPHENO <- memoise::memoise(get_pheno)

get_resdf <- function(){
    fn = "../../../residualized_expression/_m/genes_residualized_expression.csv"
    return(data.table::fread(fn) %>%
           dplyr::filter(gene_id %in% get_eGenes()$gene_id) %>%
           tibble::column_to_rownames("gene_id") %>% t %>%
           as.data.frame %>% tibble::rownames_to_column("RNum"))
}
memRES <- memoise::memoise(get_resdf)

merge_data <- function(){
    return(dplyr::inner_join(memPHENO(), memRES(), by="RNum"))
}
memDF <- memoise::memoise(merge_data)

##### MAIN
main <- function(){
    pvals <- c()
    genes <-  dplyr::select(memDF(), any_of(get_eGenes()$gene_id)) %>% colnames
    for(gene_id in genes){
        model = paste0(gene_id, " ~ Region")
        fitted = anova(lm(model, data=memDF()))
        pvals = c(pvals, fitted["Region", "Pr(>F)"])
    }
    #pval_df = data.frame("gene_id"=genes, "p_values"=pvals)
    print(sum(pvals < 0.05))
}

main()

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()
