## Author: Louise Huuki
## Edited by KJ Benjamin

library(tidyverse)

## Functions
save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

pca_norm_data <- function(){
    ## Load voom normalized data
    fname = "../../../differential_expression/_m/genes/voomSVA.RData"
    load(fname)
    ## Transpose expression
    norm_df = v$E %>% t
    ## Calculate PCA
    pca_df = prcomp(norm_df, center=TRUE)$x
    ## Convert to data frame
    norm_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        select(c(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)) %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor) %>% rename("RNum"="sample")
    return(norm_dt)
}

pca_res_data <- function(){
    ## Read in residualized data
    fname = "../../../differential_expression/_m/genes/residualized_expression.tsv"
    res_df = data.table::fread(fname) %>% column_to_rownames("V1") %>% t
    ## Calculate PCA
    pca_df = prcomp(res_df, center=TRUE)$x
    res_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        select(c(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)) %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor) %>% rename("RNum"="sample")
    return(res_dt)
}

get_pheno <- function(){
    fname = paste0("/ceph/projects/v4_phase3_paper/inputs/",
                   "phenotypes/_m/caudate_phenotypes.csv")
    df = data.table::fread(fname) %>% filter(Dx %in% c("CTL"), Age > 17)
    return(df)
}
memPHENO <- memoise::memoise(get_pheno)

get_cell_prop <- function(){
    ## Load Bisque Estimated Props
    load("../../_h/est_prop_Bisque.Rdata")
    df = est_prop_bisque$caudate$Est.prop.long %>%
        inner_join(memPHENO(), by=c("sample"="RNum"))
    return(df)
}
memPROP <- memoise::memoise(get_cell_prop)

#### Correlation with expression PCs ####
df <- memPROP() %>% mutate(RNum=sample)
## Normalized expression
norm_long <- pca_norm_data()
est_prop_norm <- inner_join(df, norm_long, by="RNum") %>%
    select(RNum, cell_type, PC, PC_values, prop)
est_prop_norm$PC <- factor(est_prop_norm$PC, levels = paste0('PC', 1:10))
## Residualized expression
res_long <- pca_res_data()
est_prop_res <- inner_join(df, res_long, by="RNum") %>%
    select(RNum, cell_type, PC, PC_values, prop)
est_prop_res$PC <- factor(est_prop_res$PC, levels = paste0('PC', 1:10))
## Calculate p-values
print("Normalized:")
prop_norm_fit <- est_prop_norm %>% group_by(cell_type, PC) %>%
    do(fitNORM = broom::tidy(lm(prop ~ PC_values, data = .))) %>%
    unnest(fitNORM) %>% filter(term == "PC_values") %>%
    mutate(p.bonf = p.adjust(p.value, "bonf"),
           p.bonf.sig = p.bonf < 0.05,
           p.bonf.cat = cut(p.bonf,
                            breaks = c(1,0.05, 0.01, 0.005, 0),
                            labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                            include.lowest = TRUE),
           p.fdr = p.adjust(p.value, "fdr"),
           log.p.bonf = -log10(p.bonf))
print(prop_norm_fit %>% count(p.bonf.cat))
print("Residualized:")
prop_res_fit <- est_prop_res %>% group_by(cell_type, PC) %>%
    do(fitRES = broom::tidy(lm(prop ~ PC_values, data = .))) %>%
    unnest(fitRES) %>% filter(term == "PC_values") %>%
    mutate(p.bonf = p.adjust(p.value, "bonf"),
           p.bonf.sig = p.bonf < 0.05,
           p.bonf.cat = cut(p.bonf,
                            breaks = c(1,0.05, 0.01, 0.005, 0),
                            labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                            include.lowest = TRUE),
           p.fdr = p.adjust(p.value, "fdr"),
           log.p.bonf = -log10(p.bonf))
print(prop_res_fit %>% count(p.bonf.cat))
## Tile plot (heatmap)
my_breaks <- c(0.05, 0.01, 0.005, 0)
tile_plot_norm <- prop_norm_fit %>%
    ggplot(aes(x = cell_type, y = PC, fill = log.p.bonf,
               label=ifelse(p.bonf.sig,format(round(log.p.bonf,1), nsmall=1), ""))) +
    geom_tile(color = "grey") +
    ggfittext::geom_fit_text(contrast = TRUE, aes(fontface="bold")) +
    viridis::scale_color_viridis(option = "magma") +
    viridis::scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                                direction=-1, limits=c(0, 40)) +
    labs(title ="p-values cell-type prop~Normalized Expression", x = 'Cell Type',
         color ="p-value Bonf\nsignificance", y="Normalized Expression (PCs)") +
    ggpubr::theme_pubr(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_img(tile_plot_norm, paste0("tilePlot_norm_prop"))
tile_plot_res <- prop_res_fit %>%
    ggplot(aes(x = cell_type, y = PC, fill = log.p.bonf,
               label=ifelse(p.bonf.sig,format(round(log.p.bonf,1), nsmall=1), ""))) +
    geom_tile(color = "grey") +
    ggfittext::geom_fit_text(contrast = TRUE, aes(fontface="bold")) +
    viridis::scale_color_viridis(option = "magma") +
    viridis::scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                                direction=-1, limits=c(0, 40)) +
    labs(title ="p-values cell-type prop~Residualized Expression", x = 'Cell Type',
         color ="p-value Bonf\nsignificance", y="Residualized Expression (PCs)") +
    ggpubr::theme_pubr(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
save_img(tile_plot_res, paste0("tilePlot_res_prop"))

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
