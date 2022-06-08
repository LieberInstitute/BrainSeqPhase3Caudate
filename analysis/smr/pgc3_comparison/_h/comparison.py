"""
This script reviews the SMR results in context with GTEx dataset.
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
def get_smr():
    fn = "../../_m/eqtl_genes.eqtl_p1e-04.gwas_p5e-08.csv"
    return pd.read_csv(fn, sep=',')\
        .merge(gene_annotation(), left_on="probeID", right_on="featureID")


@lru_cache()
def get_pgc3():
    ## Extract prioritized genes from PGC3 (FINEMAP or SMR evidence)
    gwas_fn = '/ceph/users/jbenja13/resources/gwas/pgc3/_m/'+\
        'nature_submission_11.08.2021/Supplementary Tables/'+\
        'Supplementary Table 12 - Prioritized Genes UPDATED.xlsx'
    return pd.read_excel(gwas_fn, sheet_name="Prioritised")


@lru_cache()
def merge_data():
    bs = get_smr().loc[:, ["probeID", "p_SMR", "p_HEIDI", "FDR", "Symbol", "ensemblID"]]\
                  .drop_duplicates(subset="probeID")
    pgc3 = get_pgc3()
    return pd.merge(bs, pgc3, left_on="ensemblID", right_on="Ensembl.ID",
                    suffixes=("_BS", "_PGC3"))


def main():
    df = merge_data()[(merge_data()["FDR"] < 0.05) &
                      (merge_data()["p_HEIDI"] > 0.01)].copy()
    with open("summary.log", "w") as f:
        print("There are %d shared SMR assocations (SMR < 0.05, HEIDI > 0.01)!"
              % df.shape[0], file=f)
        print(df.Symbol.values, file=f)


if __name__ == "__main__":
    main()
