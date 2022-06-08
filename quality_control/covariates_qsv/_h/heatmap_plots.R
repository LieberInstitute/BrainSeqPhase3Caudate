## Author: Louise Huuki
## Edited by KJ Benjamin

library(tidyverse)

## Functions
save_img <- function(image, fn, w=14, h=7){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

pca_norm_data <- function(){
    ## Load voom normalized data
    fname = "../../../analysis/differential_expression/_m/genes/voomSVA.RData"
    load(fname)
    ## Transpose expression
    norm_df = v$E %>% t
    ## Calculate PCA
    pca_df = prcomp(norm_df, center=TRUE)$x
    ## Convert to data frame
    norm_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        select(c(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                 PC11, PC12, PC13)) %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor) %>% rename("RNum"="sample")
    return(norm_dt)
}
memNORM <- memoise::memoise(pca_norm_data)

pca_res_data <- function(){
    ## Read in residualized data
    fname = paste0("../../../analysis/differential_expression/",
                   "_m/genes/residualized_expression.tsv")
    res_df = data.table::fread(fname) %>% column_to_rownames("V1") %>% t
    ## Calculate PCA
    pca_df = prcomp(res_df, center=TRUE)$x
    res_dt = pca_df %>% as.data.frame %>% rownames_to_column("sample") %>%
        select(c(sample, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                 PC11, PC12, PC13)) %>%
        pivot_longer(-sample, names_to="PC", values_to="PC_values") %>%
        mutate_if(is.character, as.factor) %>% rename("RNum"="sample")
    return(res_dt)
}
memRES <- memoise::memoise(pca_res_data)

get_pheno <- function(){
    fname = paste0("/ceph/projects/v4_phase3_paper/inputs/",
                   "phenotypes/_m/caudate_phenotypes_full.csv")
    df = data.table::fread(fname) %>% column_to_rownames("RNum") %>%
        filter(Dx %in% c("CTL", "SZ"), Age > 17) %>%
        mutate(across(where(is.character), as.factor)) %>%
        mutate(across(where(is.factor), as.numeric)) %>%
        mutate(across(where(is.logical), as.numeric))
    return(df)
}
memPHENO <- memoise::memoise(get_pheno)

prep_data <- function(covars){
    df = covars %>%
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
    return(df)
}
memEST <- memoise::memoise(prep_data)

prep_qsv <- function(){
    qsv_filename = paste0("/ceph/projects/v4_phase3_paper/inputs/counts/",
                          "text_files_counts/_m/caudate/qSV_caudate.csv")
    est_df <- read.csv(qsv_filename, row.names = "X") %>%
        rownames_to_column("RNum") %>%
        pivot_longer(!RNum, names_to = "qSV", values_to = "qSV_value")
    est_df$qSV <- factor(est_df$qSV, levels = paste0('PC', 1:15))
    return(est_df)
}
memQSV <- memoise::memoise(prep_qsv)

merge_data <- function(covars){
    df <- inner_join(memEST(covars), memQSV(), by="RNum")
    return(df)
}
memDF <- memoise::memoise(merge_data)

merge_expr <- function(covars, fnc){
    df <- inner_join(memEST(covars), fnc(), by="RNum")
    df$PC <- factor(df$PC, levels = paste0('PC', 1:15))
    return(df)
}
memEXPR <- memoise::memoise(merge_expr)

fit_model <- function(covars, qsv, norm){
    if(qsv){
        est_fit0 <- memDF(covars) %>% group_by(Covariate, qSV) %>%
            do(fitEST = broom::tidy(lm(Variable ~ qSV_value, data = .)))
    } else {
        if(norm){
            est_fit0 <- memEXPR(covars, memNORM) %>% group_by(Covariate, PC) %>%
                do(fitEST = broom::tidy(lm(Variable ~ PC_values, data = .)))
        } else {
            est_fit0 <- memEXPR(covars, memRES) %>% group_by(Covariate, PC) %>%
                do(fitEST = broom::tidy(lm(Variable ~ PC_values, data = .)))
        }
    }
    ## Calculate p-values
    est_fit <- est_fit0 %>% unnest(fitEST) %>% filter(term != "(Intercept)") %>%
        mutate(p.bonf = p.adjust(p.value, "bonf"), p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf, breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"), log.p.bonf = -log10(p.bonf))
    print(est_fit %>% count(p.bonf.cat))
    return(est_fit)
}

tile_plot <- function(covars, qsv=TRUE, norm=TRUE){
    ## Tile plot (heatmap)
    my_breaks <- c(0.05, 0.01, 0.005, 0)
    limits = c(0, 75)
    xlabel = "Covariate"
    if(qsv){
        ylabel = "qSV"; out = "qsv"
        tile_plot0 <- fit_model(covars, qsv, norm) %>%
            ggplot(aes(x = Covariate, y = qSV, fill = log.p.bonf,
                       label=ifelse(p.bonf.sig,format(round(log.p.bonf,1),
                                                      nsmall=1), "")))
    } else {
        tile_plot0 <- fit_model(covars, qsv, norm) %>%
            ggplot(aes(x = Covariate, y = PC, fill = log.p.bonf,
                       label=ifelse(p.bonf.sig,format(round(log.p.bonf,1),
                                                      nsmall=1), "")))
        if(norm){
            ylabel = "Normalized Expression"; out = "norm"
        } else {
            ylabel = "Residualized Expression"; out = "res"
        }
    }
    tile_plot <- tile_plot0 + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE, aes(fontface="bold")) +
        viridis::scale_color_viridis(option = "magma") +
        viridis::scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                                    direction=-1, limits=limits) +
        labs(x=xlabel, y=ylabel) + ggpubr::theme_pubr(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    save_img(tile_plot, paste0("tilePlot_",out,"_covariates"))
}

#### Correlation with expression PCs ####
covarsCont = memPHENO() %>% rownames_to_column("RNum") %>%
    select(-c("Region", "BrNum", "Protocol", "antipsychotics",
              "lifetime_antipsych", "bamFile", "IID", "SOL",
              "trimmed", "SAMPLE_ID"))
## Plotting
tile_plot(covarsCont)
tile_plot(covarsCont, FALSE, TRUE)
tile_plot(covarsCont, FALSE, FALSE)

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
