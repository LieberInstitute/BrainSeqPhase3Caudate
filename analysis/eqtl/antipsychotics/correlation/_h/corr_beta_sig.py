"""
This script examines the Spearman correlation of beta (slope) across control
only, SZ with antipsychotics medication (AP), and SZ without AP.

This uses significant SNP-gene pairs.
"""
import functools
import pandas as pd
from scipy.stats import spearmanr
from rpy2.robjects.packages import importr
from rpy2.robjects import r, globalenv, pandas2ri

@functools.lru_cache()
def load_beta(label):
    label_dict = {"control_only": "CTL", "sz_AP":"SZ_AP", "sz_noAP": "SZ_noAP"}
    fn ="/ceph/projects/v4_phase3_paper/analysis/eqtl_analysis/"+\
        "%s/genes/expression_gct/prepare_expression/" % label+\
        "annotate_outputs/_m/Brainseq_LIBD.signifpairs.txt.gz"
    return pd.read_csv(fn, sep='\t', usecols=[0,1,7])\
             .rename(columns={"slope":label_dict[label]})


@functools.lru_cache()
def prepare_data():
    ctl = load_beta("control_only")
    szAP = load_beta("sz_AP")
    noAP = load_beta("sz_noAP")
    return ctl, szAP, noAP


@functools.lru_cache()
def merge_data():
    ctl, szAP, noAP = prepare_data()
    return pd.merge(ctl, szAP, on=["gene_id", "variant_id"])\
             .merge(noAP, on=["gene_id", "variant_id"])


def corr_beta(df, label1, label2):
    return spearmanr(df[label1], df[label2])


def plotNsave_corr(df):
    #ggpubr = importr('ggpubr')
    pandas2ri.activate()
    globalenv['df'] = df
    r('''
    pp1 = ggpubr::ggscatter(df, x="CTL", y="SZ_AP", add="reg.line", size=1,
                            xlab="eQTL beta (CTL)", ylab="eQTL beta (SZ w/ AP)",
                            add.params=list(color="blue", fill="lightgray"),
                            conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                            cor.method="spearman", cor.coeff.args=list(label.sep="\n"),
                            ggtheme=ggpubr::theme_pubr(base_size=15))
    pp2 = ggpubr::ggscatter(df, x="CTL", y="SZ_noAP", add="reg.line", size=1,
                            xlab="eQTL beta (CTL)", ylab="eQTL beta (SZ w/o AP)",
                            add.params=list(color="blue", fill="lightgray"),
                            conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                            cor.method="spearman", cor.coeff.args=list(label.sep="\n"),
                            ggtheme=ggpubr::theme_pubr(base_size=15))
    pp3 = ggpubr::ggscatter(df, x="SZ_AP", y="SZ_noAP", add="reg.line", size=1,
                            xlab="eQTL beta (SZ w/ AP)", conf.int=TRUE,
                            ylab="eQTL beta (SZ w/o AP)", cor.coef=TRUE,
                            cor.coef.size=3, cor.method="spearman",
                            cor.coeff.args=list(label.sep="\n"),
                            add.params=list(color="blue", fill="lightgray"),
                            ggtheme=ggpubr::theme_pubr(base_size=15))
    figure = ggpubr::ggarrange(pp1, pp2, pp3, ncol=3, nrow=1, align='v')
    ggplot2::ggsave(plot=figure, dpi=72, width=9, height=3,
                    filename="scatter_slope_antipsychotics_eqtl_sig.pdf")
    ggplot2::ggsave(plot=figure, dpi=300, width=9, height=3,
                    filename="scatter_slope_antipsychotics_eqtl_sig.png")
    ''')


def main():
    ## Load data
    df = merge_data()
    ## Calculate rho
    with open("rho_statistics_sig.log", "w") as f:
        rho1, pval1 = corr_beta(df, "CTL", "SZ_AP")
        print("rho between CTL and SZ AP:\t\t %.3f" % rho1, file=f)
        rho2, pval2 = corr_beta(df, "CTL", "SZ_noAP")
        print("rho between CTL and SZ no AP:\t\t %.3f" % rho2, file=f)
        rho3, pval3 = corr_beta(df, "SZ_AP", "SZ_noAP")
        print("rho between SZ with and without AP:\t %.3f" % rho3, file=f)
    ## Generate figure
    plotNsave_corr(df)


if __name__ == "__main__":
    main()
