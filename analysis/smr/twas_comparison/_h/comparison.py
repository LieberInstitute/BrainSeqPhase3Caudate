"""
This script compares TWAS with SMR results.
"""
import pandas as pd
from functools import lru_cache
from scipy.stats import fisher_exact, spearmanr

@lru_cache()
def gene_annotation():
    fn = "/ceph/projects/v4_phase3_paper/inputs/counts/gene_annotation/"+\
        "_m/gene_annotation.tsv"
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def get_twas():
    """
    Load all TWAS results (heritable features).
    """
    fn = "/ceph/projects/v4_phase3_paper/analysis/twas_ea/gene_weights/"+\
        "fusion/summary_stats/_m/fusion_associations.txt"
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def get_smr():
    """
    Load SMR results (all).
    """
    fn = "../../_m/eqtl_genes.eqtl_p1e-04.gwas_p5e-08.csv"
    return pd.read_csv(fn, sep=',')\
             .merge(gene_annotation(), left_on="probeID", right_on="featureID")


@lru_cache()
def cal_enrichment():
    twas_uni = set(get_twas().FILE)
    twas_sig = set(get_twas()[(get_twas()["FDR"] < 0.05)].FILE)
    smr_uni = set(get_smr().ensemblID)
    smr_sig = set(get_smr()[(get_smr()["FDR"] < 0.05) &
                            (get_smr()["p_HEIDI"] > 0.01)].ensemblID)
    yes_twas = twas_sig; yes_smr = smr_sig
    no_twas = twas_uni - twas_sig; no_smr = smr_uni - smr_sig
    m = [[len(yes_twas.intersection(yes_smr)),len(no_twas.intersection(yes_smr))],
         [len(yes_twas.intersection(no_smr)),len(no_twas.intersection(no_smr))]]
    return fisher_exact(m)


@lru_cache()
def corr_beta():
    twas = get_twas().loc[:, ["FILE", "TWAS.Z"]]
    smr = get_smr().loc[:, ["ensemblID", "b_SMR"]]
    df = pd.merge(twas, smr, left_on="FILE", right_on="ensemblID")
    return spearmanr(df.b_SMR, df["TWAS.Z"])


def main():
    yes_twas = set(get_twas()[(get_twas()["FDR"] < 0.05)].FILE)
    yes_smr = set(get_smr()[(get_smr()["FDR"] < 0.05) &
                            (get_smr()["p_HEIDI"] > 0.01)].ensemblID)
    with open("summary.log", "w") as f:
        print("{} overlapping genes between TWAS and SMR!"\
              .format(len(yes_twas.intersection(yes_smr))), file=f)
        print("{:.1%} of TWAS also SMR associations!"\
              .format(len(yes_twas.intersection(yes_smr)) / len(yes_smr)),
              file=f)
        print("Fisher Exact Test, Odd ratio = %.2f, p-value < %.1e" %
              cal_enrichment(), file=f)
        print("Spearman correlation, rho = %.2f, p-value < %.1e" %
              corr_beta(), file=f)
    get_smr()[(get_smr()["ensemblID"].isin(yes_twas.intersection(yes_smr)))]\
        .to_csv("BrainSeq_caudate_SMR_TWAS.tsv.gz", sep='\t', index=False)


if __name__ == '__main__':
    main()
