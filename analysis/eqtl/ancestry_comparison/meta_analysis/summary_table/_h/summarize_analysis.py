"""
This script summarizes eQTL analysis.
"""
import pandas as pd
from functools import lru_cache

config = {
    "genes": "/ceph/projects/v4_phase3_paper/inputs/counts/"+\
    "text_files_counts/_m/caudate/gene_annotation.tsv",
    "transcripts": "/ceph/projects/v4_phase3_paper/inputs/counts/"+\
    "text_files_counts/_m/caudate/tx_annotation.tsv",
    "exons": "/ceph/projects/v4_phase3_paper/inputs/counts/"+\
    "text_files_counts/_m/caudate/exon_annotation.tsv",
    "junctions": "/ceph/projects/v4_phase3_paper/inputs/counts/"+\
    "text_files_counts/_m/caudate/jxn_annotation.tsv"
}

@lru_cache()
def get_mashr_eqtls(feature):
    df = pd.read_csv("../../_m/%s/lfsr_allpairs_ancestry.txt.gz" % feature,
                     sep='\t')
    ## Fix bug in compute_lfsr, where negative values should be zero
    df["AA"] = [0 if x < 0 else x for x in df.AA]
    df["EA"] = [0 if x < 0 else x for x in df.EA]
    return df[(df["AA"] < 0.05) | (df["EA"] < 0.05)]


@lru_cache()
def annotate_eqtls(feature):
    annot = pd.read_csv(config[feature], sep='\t')\
              .loc[:, ["names", "seqnames", "gencodeID"]]
    return pd.merge(get_mashr_eqtls(feature), annot,
                    left_on="gene_id", right_on="names")\
             .drop(["names"], axis=1)


@lru_cache()
def load_pgc2():
    pgc2_file = '/ceph/projects/v4_phase3_paper/inputs/sz_gwas/'+\
        'pgc2_clozuk/map_phase3/_m/libd_hg38_pgc2sz_snps_p5e_minus8.tsv'
    return pd.read_csv(pgc2_file, sep='\t', low_memory=False)


@lru_cache()
def load_pgc3():
    pgc3_file = '/ceph/projects/v4_phase3_paper/inputs/sz_gwas/'+\
        'pgc3/map_phase3/_m/libd_hg38_pgc2sz_snps_p5e_minus8.tsv'
    return pd.read_csv(pgc3_file, sep='\t', low_memory=False)


@lru_cache()
def merge_pgc2_N_eqtl(feature):
    return load_pgc2().merge(annotate_eqtls(feature), how='inner',
                             left_on='our_snp_id', right_on='variant_id',
                             suffixes=['_PGC2', '_EQTL'])


@lru_cache()
def merge_pgc3_N_eqtl(feature):
    return load_pgc3().merge(annotate_eqtls(feature), how='inner',
                             left_on='our_snp_id', right_on='variant_id',
                             suffixes=['_PGC3', '_EQTL'])


@lru_cache()
def extract_features(fnc):
    ## Extract significant eQTL using mashr
    genes = fnc("genes")
    trans = fnc("transcripts")
    exons = fnc("exons")
    juncs = fnc("junctions")
    return genes, trans, exons, juncs


@lru_cache()
def get_eQTL_result_by_ancestry(fnc):
    genes, trans, exons, juncs = extract_features(fnc)
    genes["Type"] = "Gene"
    trans["Type"] = "Transcript"
    exons["Type"] = "Exon"
    juncs["Type"] = "Junction"
    df = pd.concat([genes, trans, exons, juncs])
    df["Type"] = df.Type.astype("category").cat\
                   .reorder_categories(["Gene","Transcript","Exon","Junction"])
    return df


def output_summary(fnc, variable):
    ## Extract eQTL using mashr
    genes, trans, exons, juncs = extract_features(fnc)
    ## Total significant eQTLs
    gg = len(set(genes[variable]))
    tt = len(set(trans[variable]))
    ee = len(set(exons[variable]))
    jj = len(set(juncs[variable]))
    print("\neGene:\t\t%d\neTranscript:\t%d\neExon:\t\t%d\neJunction:\t%d" %
          (gg, tt, ee, jj))


def print_summary(fnc, label):
    genes, trans, exons, juncs = extract_features(fnc)
    if label == "mashr":
        w_mode = "w"
        statement = "Significant eQTL (lfsr < 0.05) in at least one ancestry"
    else:
        w_mode = "a"
        if label == "pgc2":
            statement = "Overlap with PGC2+CLOZUK"
        else:
            statement = "Overlap with PGC3"
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        for variable in ["effect", "gene_id", "gencodeID"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            jj = len(set(juncs[variable]))
            print("\neGene:\t\t%d\neTranscript:\t%d\neExon:\t\t%d\neJunction:\t%d" %
                  (gg, tt, ee, jj), file=f)


def main():
    # Summarize results mashr (lfsr < 0.05) in at least 1 ancestry
    print_summary(annotate_eqtls, "mashr")
    get_eQTL_result_by_ancestry(annotate_eqtls)\
        .to_csv("BrainSeq_caudate_eQTL.txt.gz", sep="\t", index=False)
    # Overlap with PGC2+CLOZUK
    print_summary(merge_pgc2_N_eqtl, "pgc2")
    get_eQTL_result_by_ancestry(merge_pgc2_N_eqtl)\
        .to_csv("BrainSeq_caudate_eQTL_pgc2.txt.gz", sep="\t", index=False)
    # Overlap with PGC3
    print_summary(merge_pgc3_N_eqtl, "pgc3")
    get_eQTL_result_by_ancestry(merge_pgc3_N_eqtl)\
        .to_csv("BrainSeq_caudate_eQTL_pgc3.txt.gz", sep="\t", index=False)


if __name__ == '__main__':
    main()
