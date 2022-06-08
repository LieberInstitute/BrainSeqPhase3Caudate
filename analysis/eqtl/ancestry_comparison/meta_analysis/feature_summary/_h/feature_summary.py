"""
This script generates a summary for mashr results by feature.
"""
import functools
import pandas as pd

@functools.lru_cache()
def get_data(feature):
    fn = "../../_m/%s/lfsr_allpairs_ancestry.txt.gz" % feature
    df = pd.read_csv(fn, sep='\t', low_memory=True)
    aa = df[(df["AA"] < 0.05) & (df["EA"] >= 0.05)].copy()
    ea = df[(df["AA"] >= 0.05) & (df["EA"] < 0.05)].copy()
    shared_df = df[(df["AA"] < 0.05) & (df["EA"] < 0.05)].copy()
    all_df = df[(df["AA"] < 0.05) | (df["EA"] < 0.05)].copy()
    return aa, ea, shared_df, all_df


def feature_summary(feature):
    aa, ea, shared_df, all_df = get_data(feature)
    print(feature.upper())
    print("There are %d AA specific SNP-feature!" % aa.shape[0])
    print("There are %d EA specific SNP-feature!" % ea.shape[0])
    print("There are {} ({:.1%}) SNP-feature shared between ancestry!\n"\
          .format(shared_df.shape[0], shared_df.shape[0] / all_df.shape[0]))


def efeature_summary(feature):
    aa, ea, shared_df, all_df = get_data(feature)
    aa = aa.groupby("gene_id").first().reset_index()
    ea = ea.groupby("gene_id").first().reset_index()
    shared_df = shared_df.groupby("gene_id").first().reset_index()
    all_df = all_df.groupby("gene_id").first().reset_index()
    print(feature.upper())
    print("There are %d AA specific eFeatures!" % aa.shape[0])
    print("There are %d EA specific eFeatures!" % ea.shape[0])
    print("There are {} ({:.1%}) eFeatures shared between ancestry!\n"\
          .format(shared_df.shape[0], shared_df.shape[0] / all_df.shape[0]))


def main():
    for feature in ["genes", "transcripts", "exons", "junctions"]:
        feature_summary(feature)
        efeature_summary(feature)


if __name__ == "__main__":
    main()
