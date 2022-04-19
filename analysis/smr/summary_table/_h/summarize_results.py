"""
This script summarizes the results from SMR analysis performed by Danny Chen.
"""

import pandas as pd
from functools import lru_cache

@lru_cache()
def get_smr_results(feature):
    fn = "../../_m/eqtl_%s.Caudate.CAUC_NC_SCZ_BIP.age13." % feature+\
        "index_p1e-04.SCZ_PGC3_p1e-04.smr_q0.05.heidi_p0.01.csv.gz"
    return pd.read_csv(fn)


@lru_cache()
def get_top_genes():
    cols = ['gencodeID', 'Symbol', 'seqnames', 'Probe_bp', 'topSNP',
            'b_SMR', 'se_SMR', 'p_SMR', 'p_HEIDI','nsnp_HEIDI', 'q_SMR']
    get_smr_results("gene").loc[:, cols].sort_values("p_SMR")\
        .groupby("gencodeID").first().reset_index().sort_values("q_SMR")\
        .head(25).to_csv("smr_top25.tsv", sep='\t', index=False)


@lru_cache()
def clean_data():
    cols = ['probeID', 'gencodeID', 'Symbol', 'seqnames', 'Probe_bp', 'topSNP',
            'A1', 'A2', 'Freq', 'b_GWAS', 'se_GWAS', 'p_GWAS', 'b_eQTL',
            'se_eQTL', 'p_eQTL', 'b_SMR', 'se_SMR', 'p_SMR', 'p_HEIDI',
            'nsnp_HEIDI', 'clm_id', 'q_SMR']
    genes = get_smr_results("gene").loc[:, cols]
    trans = get_smr_results("tx")\
        .rename(columns={"gene_name": "Symbol", "gene_id": "gencodeID"})\
        .loc[:, cols]
    exons = get_smr_results("exon").loc[:, cols]
    juncs = get_smr_results("jxn")\
        .drop(["Symbol"], axis=1)\
        .rename(columns={"newGeneID": "gencodeID", "newGeneSymbol": "Symbol"})\
        .loc[:, cols]
    return genes, trans, exons, juncs


@lru_cache()
def prep_smr():
    genes, trans, exons, juncs = clean_data()
    genes["Type"] = "Gene"
    trans["Type"] = "Transcript"
    exons["Type"] = "Exon"
    juncs["Type"] = "Junction"
    df = pd.concat([genes, trans, exons, juncs])
    df["Type"] = df.Type.astype("category").cat\
                   .reorder_categories(["Gene","Transcript","Exon","Junction"])
    return df


def print_summary():
    genes, trans, exons, juncs = clean_data()
    w_mode = "w"
    statement = "SMR results"
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        for variable in ["clm_id", "probeID", "gencodeID"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            jj = len(set(juncs[variable]))
            print("\nGene:\t\t%d\nTranscript:\t%d\nExon:\t\t%d\nJunction:\t%d\n" %
                  (gg, tt, ee, jj), file=f)


def main():
    get_top_genes()
    print_summary()
    prep_smr().to_csv("BrainSeq_caudate_smr_q0.05.txt.gz",
                      sep="\t", index=False)


if __name__ == '__main__':
    main()
