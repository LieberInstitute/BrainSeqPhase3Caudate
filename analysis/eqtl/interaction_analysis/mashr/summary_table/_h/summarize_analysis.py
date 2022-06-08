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
    return pd.read_csv("../../_m/%s/lfsr_allpairs_3tissues.txt.gz" %
                       feature, sep='\t')


@lru_cache()
def get_eqtl_tissue(feature, tissue):
    return get_mashr_eqtls(feature)[(get_mashr_eqtls(feature)[tissue] < 0.05)].copy()


@lru_cache()
def get_eqtls_all(feature):
    df = get_mashr_eqtls(feature)
    return df[(df["Caudate"] < 0.05) | (df["DLPFC"] < 0.05) |
              (df["Hippocampus"] < 0.05)].copy()


@lru_cache()
def get_eqtls(feature):
    cc = get_eqtl_tissue(feature, "Caudate")
    dd = get_eqtl_tissue(feature, "DLPFC")
    hh = get_eqtl_tissue(feature, "Hippocampus")
    return cc, dd, hh


@lru_cache()
def get_caudate_specific(feature):
    cc, dd, hh = get_eqtls(feature)
    caudate_only = list(set(cc.gene_id) - set(dd.gene_id) - set(hh.gene_id))
    return cc[(cc["gene_id"].isin(caudate_only))].copy()


@lru_cache()
def annotate_eqtls(feature):
    annot = pd.read_csv(config[feature], sep='\t')\
              .loc[:, ["names", "seqnames", "gencodeID"]]
    return pd.merge(get_caudate_specific(feature), annot,
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
    #juncs = fnc("junctions")
    return genes, trans, exons#, juncs


@lru_cache()
def get_eQTL_result_by_tissue(fnc):
    genes, trans, exons = extract_features(fnc)
    genes["Type"] = "Gene"
    trans["Type"] = "Transcript"
    exons["Type"] = "Exon"
    #juncs["Type"] = "Junction"
    df = pd.concat([genes, trans, exons])
    df["Type"] = df.Type.astype("category").cat\
                   .reorder_categories(["Gene","Transcript","Exon"])
    return df


@lru_cache()
def get_brain_specific(feature):
    cc, dd, hh = get_eqtls(feature)
    caudate_only = list(set(cc.gene_id) - set(dd.gene_id) - set(hh.gene_id))
    dlpfc_only = list(set(dd.gene_id) - set(cc.gene_id) - set(hh.gene_id))
    hippo_only = list(set(hh.gene_id) - set(cc.gene_id) - set(dd.gene_id))
    all_eqtls = list(set(cc.gene_id) | set(dd.gene_id) | set(hh.gene_id))
    return caudate_only, dlpfc_only, hippo_only, all_eqtls


def print_specificity():
    feature = "genes"
    cc, dd, hh, eqtls = get_brain_specific(feature)
    genes = cc + dd + hh
    snames = ["Caudate"] * len(cc) + ["DLPFC"] * len(dd) + ["Hippocampus"] * len(hh)
    pd.DataFrame({"gene_id": genes, "tissue": snames})\
        .to_csv("brain_region_specific_eGenes.tsv", sep='\t', index=False)
    with open("specificity.log", mode="w") as f:
        print("There are %d total eGenes across brain regions." % len(eqtls),
              file=f)
        print("There are %d (%.1f%%) caudate specific eGenes!" %
              (len(cc), len(cc)/len(eqtls)*100), file=f)
        print("There are %d (%.1f%%) DLPFC specific eGenes!" %
              (len(dd), len(dd)/len(eqtls)*100), file=f)
        print("There are %d (%.1f%%) hippocampus specific eGenes!\n" %
              (len(hh), len(hh)/len(eqtls)*100), file=f)


def print_summary(fnc, label):
    genes, trans, exons = extract_features(fnc)
    if label == "mashr":
        w_mode = "w"
        statement = "Significant eQTL (lfsr < 0.05), caudate specific"
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
            #jj = len(set(juncs[variable]))
            print("\neGene:\t\t%d\neTranscript:\t%d\neExon:\t\t%d" %
                  (gg, tt, ee), file=f)


def main():
    print_specificity()
    # Summarize results mashr (lfsr < 0.05) in at least 1 tissue
    print_summary(annotate_eqtls, "mashr")
    get_eQTL_result_by_tissue(annotate_eqtls)\
        .to_csv("BrainSeq_caudateSpecific_eQTL.txt.gz", sep="\t", index=False)
    # Overlap with PGC2+CLOZUK
    print_summary(merge_pgc2_N_eqtl, "pgc2")
    get_eQTL_result_by_tissue(merge_pgc2_N_eqtl)\
        .to_csv("BrainSeq_caudateSpecific_eQTL_pgc2.txt.gz", sep="\t",
                index=False)
    # Overlap with PGC3
    print_summary(merge_pgc3_N_eqtl, "pgc3")
    get_eQTL_result_by_tissue(merge_pgc3_N_eqtl)\
        .to_csv("BrainSeq_caudateSpecific_eQTL_pgc3.txt.gz", sep="\t",
                index=False)


if __name__ == '__main__':
    main()
