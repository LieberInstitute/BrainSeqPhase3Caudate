## Author: Louise Huuki
## Edited by KJ Benjamin

library(SummarizedExperiment)
library(DeconvoBuddies)
library(RColorBrewer)
library(sessioninfo)
library(patchwork)
library(tidyverse)
library(viridis)
library(ggpubr)
library(broom)

## Function to save plots
save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

get_pheno <- function(){
    df = data.table::fread(paste0("/ceph/projects/v4_phase3_paper/inputs/",
                                  "phenotypes/_m/merged_phenotypes.csv")) %>%
        filter(Dx %in% c("SZ", "CTL"), Age > 17)
    return(df)
}

memPHENO <- memoise::memoise(get_pheno)
## Create directory for plots
dir.create("plots")
## Load colors and plotting functions
cell_colors <- create_cell_colors(pallet = "classic")
cell_colors <- cell_colors[c("Astro", "Endo", "Micro", "Mural", "Oligo",
                             "OPC", "Tcell", "Excit", "Inhib")]
## Load Bisque Estimated Props
load("../../_h/est_prop_Bisque.Rdata", verbose = TRUE)

long_prop <- est_prop_bisque$caudate$Est.prop.long %>%
  rename(RNum = sample) %>% inner_join(memPHENO(), by="RNum")

## Boxplots
group_boxPlot <- long_prop %>%
  ggplot(aes(x = cell_type, y = prop, fill = cell_type)) +
  geom_boxplot() + labs(x = "Cell Type", y = "Proportion") +
  theme_pubr(base_size = 15) +
  scale_fill_manual(values = cell_colors, guide = "none")

save_img(group_boxPlot, "plots/cellType_boxplots", h=6)

#### composition barplot ####
comp_barplot <- plot_composition_bar(long_prop, x_col="Dx") +
    theme_pubr(base_size = 15) +
    scale_fill_manual(values = cell_colors)

save_img(comp_barplot, "plots/composition_barplot")

#### Cor with qSV ####
qsv_filename = paste0("/ceph/projects/v4_phase3_paper/inputs/counts/",
                      "text_files_counts/_m/caudate/qSV_caudate.csv")
qSV_mat <- read.csv(qsv_filename, row.names = "X")
dim(qSV_mat)

identical(rownames(est_prop_bisque$caudate$bulk.props), rownames(qSV_mat))

qSV_long <- qSV_mat %>% rownames_to_column("RNum") %>%
  pivot_longer(!RNum, names_to = "qSV", values_to = "qSV_value")

## Bind with qSV table
est_prop_qsv <- left_join(long_prop, qSV_long, by = "RNum")
est_prop_qsv$qSV <- factor(est_prop_qsv$qSV,
                           levels = paste0('PC', 1:13))
levels(est_prop_qsv$qSV)

#### Calculate p-values ####
prop_qSV_fit <- est_prop_qsv %>% group_by(cell_type, qSV) %>%
  do(fitQSV = tidy(lm(prop ~ qSV_value, data = .))) %>%
  unnest(fitQSV) %>% filter(term == "qSV_value") %>%
  mutate(p.bonf = p.adjust(p.value, "bonf"),
         p.bonf.sig = p.bonf < 0.05,
         p.bonf.cat = cut(p.bonf,
                          breaks = c(1,0.05, 0.01, 0.005, 0),
                          labels = c("<= 0.005","<= 0.01", "<= 0.05", "> 0.05"),
                          include.lowest = TRUE),
         p.fdr = p.adjust(p.value, "fdr"),
         log.p.bonf = -log10(p.bonf))

print(prop_qSV_fit %>% count(p.bonf.cat))
# p.bonf.cat     n
# <fct>      <int>
# 1 <= 0.005      38
# 2 <= 0.01        1
# 3 <= 0.05        5
# 4 > 0.05        73

#### Tile plots ####
my_breaks <- c(0.05, 0.01, 0.005, 0)

# sig_colors <- c(rev(viridis_pal(option = "magma")(3)),NA)
# names(sig_colors) <- levels(prop_qSV_fit$p.bonf.cat)

tile_plot_val <- prop_qSV_fit %>%
    ggplot(aes(x = cell_type, y = qSV, fill = log.p.bonf,
               label=ifelse(p.bonf.sig, format(round(log.p.bonf,1), nsmall=1), ""))) +
    geom_tile(color = "grey") +
    ggfittext::geom_fit_text(contrast = TRUE, aes(fontface="bold")) +
    scale_color_viridis(option = "magma") +
    scale_fill_viridis(name="-log10(p-value Bonf)", option="magma",
                       direction=-1) +
    labs(title ="p-values cell-type prop~qSV", x = 'Cell Type',
         color ="p-value Bonf\nsignificance") +
    theme_pubr(base_size = 15) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

save_img(tile_plot_val, "plots/qSV_prop_fit_tileVal")

#### Create scatter plots ####
sig_colors2 <- c("#440154","#31688E", "#35B779","black")
names(sig_colors2) <- levels(prop_qSV_fit$p.bonf.cat)

est_prop_qsv_fit <- left_join(est_prop_qsv, prop_qSV_fit)

scatter_plot_cat <- est_prop_qsv_fit %>%
  ggplot(aes(x = qSV_value, y = prop, color = p.bonf.cat))+
  geom_point(size = .4, alpha = .5) +
  facet_grid(cell_type~qSV, scales = "free") +
  scale_color_manual(values = sig_colors2) +
  theme_pubr(base_size = 15)+
  theme(legend.text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(size=5)))

save_img(scatter_plot_cat, "plots/qSV_cellType_scatter_cat", 26, 18)

scatter_plot_val <- est_prop_qsv_fit %>%
  ggplot(aes(x = qSV_value, y = prop, color = log.p.bonf))+
  geom_point(size = .4, alpha = .6) +
  facet_grid(cell_type~qSV, scales = "free")+
  scale_color_viridis(name = "-log10(p-value Bonf)", option = "magma", direction = -1) +
  theme_pubr(base_size = 15)+
  theme(legend.text = element_text(size = 15)) +
  guides(color = guide_legend(override.aes = list(size=5)))

save_img(scatter_plot_val, "plots/qSV_cellType_scatter_val", 26, 18)

## sgejobs::job_single('deconvo_plots', create_shell = TRUE, queue= 'bluejay',
## memory = '10G', command = "Rscript deconvo_plots.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
