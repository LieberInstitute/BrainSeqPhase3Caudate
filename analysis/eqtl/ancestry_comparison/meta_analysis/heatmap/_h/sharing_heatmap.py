"""
This script generates heatmap for eQTL sharing with mashr outputs.
"""
import functools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from rpy2.robjects import r, pandas2ri, globalenv

@functools.lru_cache()
def load_data(feature, factor):
    pandas2ri.activate()
    globalenv["feature"] = feature
    globalenv["factor"] = factor
    mat = r('''
    fn = paste0("../../_m/",feature,"/mashr_meta_results.RData")
    load(fn); mashr::get_pairwise_sharing(m2, factor=factor)
    ''')
    return pd.DataFrame(np.array(mat), index=["AA", "EA"], columns=["AA", "EA"])


def plot_heatmap(feature, factor, label):
    df = load_data(feature, factor)
    sns.set(font_scale=2)
    grid_kws = {"height_ratios": (.9, .05), "hspace": .3}
    f, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws, figsize=(10,10))
    chart = sns.heatmap(df, ax=ax, vmin=0, vmax=1, annot=True, cbar_ax=cbar_ax,
                        cbar_kws={"orientation": "horizontal"})
    chart.set_yticklabels(chart.get_yticklabels(), fontweight="bold")
    chart.set_xticklabels(chart.get_xticklabels(), fontweight="bold")
    sns_plot = chart.get_figure()
    sns_plot.savefig("eQTL_sharing_heatmap_%s_%s.pdf" % (label, feature))
    sns_plot.savefig("eQTL_sharing_heatmap_%s_%s.png" % (label, feature))


def main():
    config = {0:"signOnly", 0.5:"general", 0.95:"within_95", 0.99:"within_99"}
    for feature in ["genes", "transcripts", "exons", "junctions"]:
        for factor in [0, 0.5, 0.95, 0.99]:
            label = config[factor]
            #print(feature, label)
            plot_heatmap(feature, factor, label)


if __name__ == "__main__":
    main()
