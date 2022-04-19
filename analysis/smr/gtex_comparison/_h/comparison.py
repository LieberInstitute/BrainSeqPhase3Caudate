"""
This script reviews the SMR results in context with GTEx dataset.
"""
import pandas as pd
from functools import lru_cache
from scipy.stats import spearmanr

@lru_cache()
def get_smr():
    fn = "../../_m/eqtl_gene.Caudate.CAUC_NC_SCZ_BIP.age13.index_p1e-04"+\
        ".SCZ_PGC3_p1e-04.smr_q0.05.heidi_p0.01.csv.gz"
    return pd.read_csv(fn, sep=',')


@lru_cache()
def get_gtex():
    fn = "../_h/eqtl_gene.age13.Brain_Caudate_basal_ganglia.kb_500.r2_0.1"+\
        ".index_p1e-04.clumped.smr.csv.gz"
    return pd.read_csv(fn, sep=',')


@lru_cache()
def merge_data():
    bs = get_smr().loc[:, ["probeID", "b_SMR", "p_SMR", "p_HEIDI", "q_SMR"]]\
                  .drop_duplicates(subset="probeID")
    gtex = get_gtex().loc[:, ["probeID", "b_SMR", "p_SMR", "p_HEIDI"]]\
                     .drop_duplicates(subset="probeID")
    return pd.merge(bs, gtex, on="probeID", suffixes=("_BS", "_GTEx"))


def corr_smr():
    df = merge_data()
    return spearmanr(df.b_SMR_BS, df.b_SMR_GTEx)


def main():
    gtex = get_gtex()[(get_gtex()["p_HEIDI"] > 0.01) &
                      (get_gtex()["p_SMR"] < 0.05)].copy()
    bs = get_smr()[(get_smr()["p_HEIDI"] > 0.01) &
                   (get_smr()["q_SMR"] < 0.05)].copy()
    with open("summary.log", "w") as f:
        shared = len(set(gtex.probeID) & set(bs.probeID))
        print("There are %d shared SMR assocations (SMR < 0.05, HEIDI > 0.01)!"
              % shared, file=f)
        print("Spearman corr, BS and GTEx caudate:\nRho > %.2f, p-value < %.1e"%
              (corr_smr()), file=f)


if __name__ == "__main__":
    main()
