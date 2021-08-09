import functools
import numpy as np
import pandas as pd
from glob import iglob
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

@functools.lru_cache()
def load_pgc2_df():
    """
    Load PGC2+CLOZUK GWAS
    """
    pgc2_file = '../../../../../../inputs/sz_gwas/pgc2_clozuk/map_phase3/'+\
        '_m/libd_hg38_pgc2sz_snps.tsv'
    return pd.read_csv(pgc2_file, sep='\t', low_memory=False, index_col=0)


def load_twas_assocations(pattern, MULTI=True):
    ## Load fusion association without MHC
    li = []
    for filename in iglob(pattern):
        li.append(pd.read_csv(filename, sep='\t'))
    df = pd.concat(li, axis=0, ignore_index=True).drop(['PANEL'], axis=1)
    df['FILE'] = df.FILE.str.replace('../../_m/WEIGHTS/twas_gene_expression.',
                                     '', regex=True)\
                            .str.replace('\\..*', '', regex=True)
    if MULTI:
        ## Drop NAs and save
        df = df[~(df['TWAS.P'].isna())].copy()
        ## Add multiple testing correction
        pv1 = multipletests(df.loc[:, 'TWAS.P'], method='fdr_bh')
        pv2 = multipletests(df.loc[:, 'TWAS.P'], method='bonferroni')
        df['FDR'] = pv1[1]
        df['Bonferroni'] = pv2[1]
        ## Print results
        print("There are %d features with significant p-values." %
              np.sum(df['TWAS.P'] <= 0.05))
        print("There are %d features with significant FDR." %
              np.sum(df['FDR'] <= 0.05))
        print("There are %d features with significant Bonferroni." %
              np.sum(df['Bonferroni'] <= 0.05))
    else:
        print("There are %d features with joint analysis TWAS associations." %
              np.sum(df['JOINT.P'] <= 0.05))
    return df


@functools.lru_cache()
def merge_data(df):
    ## Merge with PGC2 dataframe
    return pd.merge(df, load_pgc2_df(), left_on='BEST.GWAS.ID',
                    right_on='our_snp_id', suffixes=['_TWAS', '_PGC2'])


def cal_enrichment(df, pval_label, MERGE=True):
    if MERGE:
        dft = merge_data(df)
    else:
        dft = df
    ## Calculate enrichment
    table =  [[np.sum((dft['P']<5e-8) & ((dft[pval_label]<.05))),
               np.sum((dft['P']<5e-8) & ((dft[pval_label]>=.05)))],
              [np.sum((dft['P']>=5e-8) & ((dft[pval_label]<.05))),
               np.sum((dft['P']>=5e-8) & ((dft[pval_label]>=.05)))]]
    print(table)
    print(fisher_exact(table))


def twas_overlaps(df):
    twas = df[(df['TWAS.P'] <= 0.05)].copy()
    snps_twas = set(load_pgc2_df()['our_snp_id']) & set(twas['BEST.GWAS.ID'])
    snps_pgc2 = set(load_pgc2_df()[(load_pgc2_df()['P']<5e-8)].loc[:,'our_snp_id']) & \
        set(twas['BEST.GWAS.ID'])
    overlap = len(snps_twas)
    overlap_sig = len(snps_pgc2)
    print('There are %0.2f%% overlap between PGC2+CLOZUK and TWAS.' %
          (overlap/len(twas.loc[:, 'BEST.GWAS.ID'].unique()) * 100))
    print('There are %0.2f%% overlap between significant PGC2+CLOZUK and TWAS.' %
          (overlap_sig/len(twas.loc[:, 'BEST.GWAS.ID'].unique())* 100))
    print("There are %d novel unique SNPs associations with SZ." %
          (overlap - overlap_sig))


def main():
    print("No MHC:")
    df_noMHC = load_twas_assocations("../../_m/*.dat")
    df_noMHC.to_csv("fusion_associations_noMHC.txt", sep="\t", index=False)
    for pval_label in ["TWAS.P", "FDR", "Bonferroni"]:
        print(pval_label)
        cal_enrichment(df_noMHC, pval_label)
    twas_overlaps(df_noMHC)
    print("All associations (including MHC):")
    df = load_twas_assocations("../../_m/*.dat*")
    df.to_csv("fusion_associations.txt", sep="\t", index=False)
    for pval_label in ["TWAS.P", "FDR", "Bonferroni"]:
        print(pval_label)
        cal_enrichment(df, pval_label)
    twas_overlaps(df)
    print("Joint associations:")
    df2 = load_twas_assocations("../../_m/sig_analysis/*included.dat", False)
    df2.to_csv("fusion_twas_joint_assoc.txt", sep="\t", index=False)
    print("JOINT.P")
    dft = df2.drop(['TWAS.Z', 'TWAS.P'], axis=1)\
             .merge(merge_data(df), on=['FILE', 'ID'], how='right').fillna(1)
    cal_enrichment(dft, "JOINT.P", False)


if __name__ == '__main__':
    main()
