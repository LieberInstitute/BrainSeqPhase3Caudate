#!/usr/bin/python

import pandas as pd
from gtfparse import read_gtf
from pybiomart import Dataset


def run_biomart():
    dataset = Dataset(name='hsapiens_gene_ensembl',
                      host='http://www.ensembl.org')
    query = dataset.query(attributes=['ensembl_gene_id',
                                      'external_gene_name'],
                          use_attr_names=True)
    biomart = query.groupby('ensembl_gene_id').first().reset_index()
    return biomart


def expr_pos():
    gtf_file = "/ceph/genome/human/gencode25/gtf.CHR/_m/gencode.v25.annotation.gtf"
    gtf0 = read_gtf(gtf_file)
    gtf = gtf0[gtf0["feature"] == "gene"].set_index("gene_id")
    return gtf[["seqname", "start", "end"]].reset_index()


def main():
    biomart = run_biomart()
    pos = expr_pos()
    pos['ensembl_gene_id'] = pos.gene_id.str.replace('\\.\d+', '', regex=True)
    df0 = pd.merge(pos, biomart, on='ensembl_gene_id')
    df0['chr'] = df0.seqname.str.replace('chr', '', regex=True)
    df = df0[['chr', 'start', 'end', 'external_gene_name']]
    df.to_csv('glist-hg38', sep=' ', index=False, header=False)


if __name__ == '__main__':
    main()
