"""
This script is used to prepare the FastQTL data for mashr.
"""
import argparse
import pandas as pd


def load_eqtl(filename):
    df = pd.read_csv(filename, sep='\t', nrows=100)
    return pd.read_csv(filename, sep='\t', dtype=df.dtypes.to_dict())


def extract_eqtls(feature):
    ## Load eQTLs for mashr
    ### AA self-reported
    aa_file = "/ceph/projects/v4_phase3_paper/analysis/eqtl_analysis/"+\
        "aa_only/%s/expression_gct/prepare_expression/" % feature +\
        "fastqtl_nominal/_m/Brainseq_LIBD.allpairs.txt.gz"
    aa_eqtl = load_eqtl(aa_file)
    ### EA self-reported
    ea_file = "/ceph/projects/v4_phase3_paper/analysis/eqtl_analysis/"+\
        "ea_only/%s/expression_gct/prepare_expression/" % feature +\
        "fastqtl_nominal/_m/Brainseq_LIBD.allpairs.txt.gz"
    ea_eqtl = load_eqtl(ea_file)
    fn = "/ceph/projects/v4_phase3_paper/analysis/eqtl_analysis/all/"+\
        "%s/expression_gct/prepare_expression/fastqtl_nominal/" % feature +\
        "_m/Brainseq_LIBD.allpairs.txt.gz"
    eqtl = load_eqtl(fn)
    return aa_eqtl, ea_eqtl, eqtl


def extract_dataframe(aa_eqtl, ea_eqtl, eqtl, variable, label, feature):
    ## AA self-reported
    df1 = aa_eqtl.loc[:, ["gene_id","variant_id",variable]]\
                 .rename(columns={variable: "AA"})
    ## EA self-reported
    df2 = ea_eqtl.loc[:, ["gene_id","variant_id",variable]]\
                 .rename(columns={variable: "EA"})
    df3 = eqtl.loc[:, ["gene_id","variant_id",variable]]\
              .rename(columns={variable: "ALL"})
    ## Merge data
    df = df1.merge(df2, on=["gene_id", "variant_id"])\
            .merge(df3, on=["gene_id", "variant_id"])
    df.to_csv("%s/%s_fastqtl_ancestry.tsv" % (feature, label),
              sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str)
    args=parser.parse_args()
    ## Main
    aa, ea, eqtl = extract_eqtls(args.feature)
    extract_dataframe(aa, ea, eqtl, "pval_nominal", "pvalue", args.feature)
    extract_dataframe(aa, ea, eqtl, "slope", "bhat", args.feature)
    extract_dataframe(aa, ea, eqtl, "slope_se", "shat", args.feature)


if __name__=='__main__':
    main()
