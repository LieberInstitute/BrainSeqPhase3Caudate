"""
This script is used to prepare GCT files for the GTEx eQTL pipeline.
Specifically for converting BrainSeq genotypes using TOPMed imputation
as well as the TPM and count files provided by BrainSeq.

Original Author: Apuã Paquola
Edited by Kynon J Benjamin for publication and modularity.

- Inputs:
  * Counts (TXT)
  * TPM (CSV)
  * Phenotype (CSV)
  * Genotype FAM information
  * Brain region of interest
- Outputs:
  * GCT files for counts and TPM for selected genes and samples
  * A lookup table of sample_ids and brain_ids
  * A list of chromosomes to use
"""
import argparse
import pandas as pd


def to_gct(filename, df):
    description_df = pd.DataFrame({'Description': df.index.values}, index=df.index)
    dfo = pd.concat([description_df, df], axis=1)
    dfo.index.name = 'Names'
    with open(filename, "wt") as out:
        print("#1.2", file=out)
        print(df.shape[0], df.shape[1], sep="\t", file=out)
        dfo.to_csv(out, sep="\t")


def tissue_map(tissue):
    return {"DLPFC": "dlpfc", "HIPPO": "hippocampus",
            "Caudate": "caudate"}[tissue]


def load_data(args):
    fam_df = pd.read_csv(args.fam_file, sep="\t", header=None,
                         names=["BrNum", "ID", "V2", "V3", "V4", "V5"])
    pheno_df = pd.read_csv(args.pheno_file)
    pheno_df = pheno_df[(pheno_df["Region"]==tissue_map(args.tissue))].copy()
    pheno_df["ids"] = pheno_df.RNum
    pheno_df.set_index("ids", inplace=True)
    tpm_df = pd.read_csv(args.tpm_file, index_col=0)
    counts_df = pd.read_csv(args.counts_file, sep="\t", index_col=0)
    return fam_df, pheno_df, tpm_df, counts_df


def select_idv(pheno_df, fam_df):
    samples = list(set(pheno_df.loc[counts_df.columns,:].BrNum)\
                   .intersection(set(fam_df.BrNum)))
    new_fam = fam_df.loc[:, (fam_df["BrNum"].isin(samples))]\
                    .drop_duplicates(subset="BrNum")
    new_fam.to_csv("keepFam.txt", sep='\t', index=False, header=False)
    return pheno_df.loc[:, ["RNum", "BrNum"]]\
                   .reset_index().set_index("BrNum")\
                   .loc[new_fam.BrNum].reset_index().set_index("ids")


def select_genes(counts_df, tpm_df):
    return list(set(counts_df.index).intersection(set(tpm_df.index)))


def generate_chrs():
    pd.DataFrame({'chr':['chr'+xx for xx in [str(x) for x in range(1,23)]+['X']]})\
      .to_csv('vcf_chr_list.txt', header=False, index=None)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pheno_file', type=str)
    parser.add_argument('--fam_file', type=str)
    parser.add_argument('--tpm_file', type=str)
    parser.add_argument('--counts_file', type=str)
    parser.add_argument('--tissue', type=str)
    args=parser.parse_args()

    generate_chrs()
    fam_df, pheno_df, tpm_df, counts_df = load_data(args)
    new_pheno = select_idv(pheno_df, fam_df)
    genes = select_genes(counts_df, tpm_df)
    to_gct("counts.gct", counts_df.loc[genes,new_pheno.index])
    to_gct("tpm.gct", tpm_df.loc[genes,new_pheno.index])
    new_pheno.loc[:, ["RNum", "BrNum"]]\
             .to_csv("sample_id_to_brnum.tsv", sep="\t", index=False)


if __name__ == '__main__':
    main()
