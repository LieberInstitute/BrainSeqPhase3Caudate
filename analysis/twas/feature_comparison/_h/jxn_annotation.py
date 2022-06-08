#!/usr/bin/python

__author__ = "KJ Benjamin"

import subprocess
import pandas as pd
from functools import lru_cache
from rpy2.robjects import r, pandas2ri, globalenv

def map_features(feature):
    return {"genes": "gene", "transcripts": "transcript",
            "exons": "exon", "junctions": "junction"}[feature]

@lru_cache()
def gene_expression(feature):
    pandas2ri.activate()
    globalenv['feature'] = feature
    r('''
    vDat = paste0("../../../differential_expression/_m/",feature,"/voomSVA.RData")
    load(vDat)
    expr = data.frame(v$E)
    expr['Geneid'] = rownames(v$E)
    ''')
    pheno_file = '/ceph/projects/v4_phase3_paper/inputs/phenotypes/_m/caudate_phenotypes.csv'
    pheno0 = pd.read_csv(pheno_file, index_col=0)
    pheno = pheno0[(pheno0['Dx'].isin(['CTL', 'SZ'])) &
                   (pheno0['Age'] > 17) & (pheno0["Race"] == "EA")].copy()
    expr0 = r['expr'].set_index('Geneid').transpose()
    expr = pheno.merge(expr0, left_index=True, right_index=True)\
                .reset_index().set_index('BrNum').drop(['index'], axis=1)\
                .drop(pheno.drop('BrNum', axis=1).columns, axis=1).transpose()
    fn = "/ceph/projects/v4_phase3_paper/inputs/counts/text_files_counts/_m/caudate/%s.bed" % map_features(feature)
    annot = pd.read_csv(fn, sep='\t', index_col=0)
    annot['chr'] = annot.seqnames.str.replace('chr', '')
    annot = annot[['gene_id', 'chr', 'start']]
    expr = pd.merge(annot, expr, left_index=True, right_index=True)
    if feature == 'junctions':
        expr = expr.reset_index().drop('index', axis=1)
        expr['JxnID'] = 'j'+expr.index.astype(str)
        expr.set_index('JxnID', inplace=True)
    return expr


def generate_annotation():
    ## Gene annotation
    jxn_df = gene_expression("junctions").loc[:, ["gene_id"]].reset_index()
    fn = "/ceph/projects/v4_phase3_paper/inputs/counts/text_files_counts/_m/" +\
        "caudate/jxn_annotation.tsv"
    annot = jxn_df.merge(pd.read_csv(fn, sep='\t', index_col=0),
                        left_on="gene_id", right_index=True)
    annot["ensemblID"] = annot.gencodeID.str.replace("\\..*", "", regex=True)
    return annot


def main():
    generate_annotation()\
        .to_csv('jxn_annotation.tsv', sep='\t', index=False, header=True)


if __name__ == '__main__':
    main()
