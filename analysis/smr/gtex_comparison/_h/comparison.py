"""
This script reviews the SMR results in context with GTEx dataset.
"""
import pandas as pd
from functools import lru_cache
from scipy.stats import fisher_exact, spearmanr

@lru_cache()
def get_smr():
    fn = "../../_m/eqtl_genes.eqtl_p1e-04.gwas_p5e-08.csv"
    return pd.read_csv(fn, sep=',')


@lru_cache()
def get_gtex():
    fn = "../_h/eqtl_gene.age13.Brain_Caudate_basal_ganglia.kb_500.r2_0.1"+\
        ".index_p1e-04.clumped.smr.csv.gz"
    return pd.read_csv(fn, sep=',')


@lru_cache()
def merge_data():
    bs = get_smr().loc[:, ["probeID", "b_SMR", "p_SMR", "p_HEIDI", "FDR"]]\
                  .drop_duplicates(subset="probeID")
    gtex = get_gtex().loc[:, ["probeID", "b_SMR", "p_SMR", "p_HEIDI"]]\
                     .drop_duplicates(subset="probeID")
    return pd.merge(bs, gtex, on="probeID", suffixes=("_BS", "_GTEx"))


@lru_cache()
def cal_enrichment():
    gtex_uni = set(get_gtex().probeID)
    gtex_sig = set(get_gtex()[(get_gtex()["p_HEIDI"] > 0.01) &
                              (get_gtex()["p_SMR"] < 0.05)].probeID)
    smr_uni = set(get_smr().probeID)
    smr_sig = set(get_smr()[(get_smr()["FDR"] < 0.05) &
                            (get_smr()["p_HEIDI"] > 0.01)].probeID)
    yes_gtex = gtex_sig; yes_smr = smr_sig
    no_gtex = gtex_uni - gtex_sig; no_smr = smr_uni - smr_sig
    m = [[len(yes_gtex.intersection(yes_smr)),len(no_gtex.intersection(yes_smr))],
         [len(yes_gtex.intersection(no_smr)),len(no_gtex.intersection(no_smr))]]
    return fisher_exact(m)


def corr_smr():
    df = merge_data()
    return spearmanr(df.b_SMR_BS, df.b_SMR_GTEx)


def main():
    gtex = get_gtex()[(get_gtex()["p_HEIDI"] > 0.01) &
                      (get_gtex()["p_SMR"] < 0.05)].copy()
    bs = get_smr()[(get_smr()["p_HEIDI"] > 0.01) &
                   (get_smr()["FDR"] < 0.05)].copy()
    with open("summary.log", "w") as f:
        shared = len(set(gtex.probeID) & set(bs.probeID))
        print("There are %d shared SMR assocations (SMR < 0.05, HEIDI > 0.01)!"
              % shared, file=f)
        print("Spearman corr, BS and GTEx caudate:\nRho > %.2f, p-value < %.1e"%
              (corr_smr()), file=f)
        print("Fisher Exact Test, Odd ratio = %.2f, p-value < %.1e" %
              cal_enrichment(), file=f)


if __name__ == "__main__":
    main()
