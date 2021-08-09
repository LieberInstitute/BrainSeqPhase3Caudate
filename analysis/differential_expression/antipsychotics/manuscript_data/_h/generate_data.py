import pandas as pd


def feature_map(feature):
    return {"genes": "Gene", "transcripts": "Transcript",
            "exons": "Exon", "junctions": "Junction"}[feature]


def get_de(feature):
    cols = ['Feature', 'gencodeID', 'ensemblID', 'Symbol', 'logFC',
            'AveExpr', 't', 'P.Value', 'adj.P.Val', 'Analysis']
    noAP = pd.read_csv("../../_m/%s/diffExpr_sz_noAPVctl_full.txt" % feature,
                       sep='\t', index_col=0)
    noAP["Feature"] = noAP.index
    noAP["Analysis"] = "SZ (no AP)vs CTL"
    AP = pd.read_csv("../../_m/%s/diffExpr_sz_APVctl_full.txt" % feature,
                     sep='\t', index_col=0)
    AP["Feature"] = AP.index
    AP["Analysis"] = "SZ (AP) vs CTL"
    if feature == "transcripts":
        rename_dict = {"gene_id": "gencodeID", "gene_name": "Symbol"}
        df = pd.concat([AP, noAP], axis=0).rename(columns=rename_dict)
        df["ensemblID"] = df.gencodeID.str.replace("\\..*", "", regex=True)
        df = df.loc[:, cols]
    elif feature == "junctions":
        rename_dict = {"newGeneID": "gencodeID", "newGeneSymbol": "Symbol"}
        df = pd.concat([AP, noAP], axis=0)\
               .drop(["Symbol"], axis=1)\
               .rename(columns=rename_dict).loc[:, cols]
    else:
        df = pd.concat([AP, noAP], axis=0).loc[:, cols]
    df["Type"] = feature_map(feature)
    return df


def main():
    lt = []
    for feature in ["genes", "transcripts", "exons", "junctions"]:
        lt.append(get_de(feature))
    pd.concat(lt).to_csv("BrainSeq_Phase3_Caudate_DiffExpr_antipsychotics_all.txt.gz",
                         sep='\t', index=False)


if __name__ == '__main__':
    main()
