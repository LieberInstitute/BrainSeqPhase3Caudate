"""
This script summarizes the results from SMR analysis performed by Danny Chen.
"""
import pandas as pd
from functools import lru_cache

@lru_cache()
def get_gene_annotation(feature):
    fn = "/ceph/projects/v4_phase3_paper/inputs/counts/gene_annotation/"+\
        "_m/%s_annotation.tsv" % feature
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def get_smr_results(feature):
    feat_map = {"gene": "genes", "tx": "transcripts",
                "exon": "exons", "jxn": "junctions"}
    fn = "../../_m/eqtl_%s.eqtl_p1e-04.gwas_p5e-08.csv" % feat_map[feature]
    return pd.read_csv(fn).merge(get_gene_annotation(feature),
                                 left_on="probeID", right_on="featureID")


@lru_cache()
def get_top_genes():
    cols = ['gencodeID', 'Symbol', 'ProbeChr', 'Probe_bp', 'topSNP', 'b_SMR',
            'se_SMR', 'p_SMR', 'p_SMR_multi', 'p_HEIDI','nsnp_HEIDI', 'FDR']
    get_smr_results("gene").loc[:, cols].sort_values("p_SMR")\
        .groupby("gencodeID").first().reset_index().sort_values("FDR")\
        .head(25).to_csv("smr_top25.tsv", sep='\t', index=False)


@lru_cache()
def clean_data():
    cols = ['probeID', 'gencodeID', 'Symbol', 'ProbeChr', 'Probe_bp', 'topSNP',
            'A1', 'A2', 'Freq', 'b_GWAS', 'se_GWAS', 'p_GWAS', 'b_eQTL',
            'se_eQTL', 'p_eQTL', 'b_SMR', 'se_SMR', 'p_SMR', 'p_SMR_multi',
            'p_HEIDI', 'nsnp_HEIDI', 'FDR']
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
def get_sig_smr(fdr=0.05):
    genes, trans, exons, juncs = clean_data()
    genes = genes[(genes["FDR"] < fdr)].copy()
    trans = trans[(trans["FDR"] < fdr)].copy()
    exons = exons[(exons["FDR"] < fdr)].copy()
    juncs = juncs[(juncs["FDR"] < fdr)].copy()
    return genes, trans, exons, juncs


@lru_cache()
def prep_smr():
    genes, trans, exons, juncs = get_sig_smr()
    genes["Type"] = "Gene"
    trans["Type"] = "Transcript"
    exons["Type"] = "Exon"
    juncs["Type"] = "Junction"
    df = pd.concat([genes, trans, exons, juncs])
    df["Type"] = df.Type.astype("category").cat\
                   .reorder_categories(["Gene","Transcript","Exon","Junction"])
    return df


def print_summary():
    genes, trans, exons, juncs = get_sig_smr()
    w_mode = "w"
    statement = "SMR results"
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        for variable in ["probeID", "gencodeID"]:
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
    prep_smr().to_csv("BrainSeq_caudate_smr_fdr0.05.txt.gz",
                      sep="\t", index=False)


if __name__ == '__main__':
    main()
