"""
This script reviews the SMR results in context with GTEx dataset.
"""
import pandas as pd
from functools import lru_cache
from scipy.stats import fisher_exact, spearmanr

@lru_cache()
def get_twas():
    fn = "../../feature_comparison/manuscript_supp_data/_m/"+\
        "BrainSeq_Phase3_Caudate_TWAS_associations_allFeatures.txt.gz"
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def get_pgc3():
    ## Extract prioritized genes from PGC3 (FINEMAP or SMR evidence)
    gwas_fn = '/ceph/users/jbenja13/resources/gwas/pgc3/_m/'+\
        'nature_submission_11.08.2021/Supplementary Tables/'+\
        'Supplementary Table 12 - Prioritized Genes UPDATED.xlsx'
    return pd.read_excel(gwas_fn, sheet_name="Prioritised")


@lru_cache()
def merge_data():
    bs = get_twas()[(get_twas()["Type"] == "Gene") &
                    (get_twas()["FDR"] < 0.05)].copy()
    pgc3 = get_pgc3()
    return pd.merge(bs, pgc3, left_on="ensemblID", right_on="Ensembl.ID",
                    suffixes=("_BS", "_PGC3"))


def main():
    df = merge_data()
    df.to_csv("BrainSeq_Phase3_Caudate_TWAS_PGC3_overlap.txt",
              sep='\t', index=False)
    with open("summary.log", "w") as f:
        print("There are %d shared TWAS assocations (FDR < 0.05)!"
              % df.shape[0], file=f)
        print(df.ID.values, file=f)


if __name__ == "__main__":
    main()
