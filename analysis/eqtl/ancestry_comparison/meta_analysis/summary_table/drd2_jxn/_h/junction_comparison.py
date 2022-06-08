"""
This script examines junction 5-7 for DRD2 between AA and EA.
"""
import pandas as pd
from functools import lru_cache

@lru_cache()
def load_data():
    cols = ["effect", "gene_id", "variant_id",
            "AA", "EA", "chr", "gencodeID", "type"]
    return pd.read_csv("junctions.txt", sep='\t', header=None, names=cols)


@lru_cache()
def extract_jxn57():
    return load_data()[(load_data()["gene_id"]=="chr11:113412884-113415420(-)")]


def main():
    aa = extract_jxn57()[(extract_jxn57()['AA'] < 0.05)].copy()
    ea = extract_jxn57()[(extract_jxn57()['EA'] < 0.05)].copy()
    shared = list(set(aa.variant_id) & set(ea.variant_id))
    with open("junction_5_7_ancestry_variants.log", "w") as f:
        print("There are a total of %d eQTL for AA" % aa.shape[0], file=f)
        print("There are a total of %d eQTL for EA" % ea.shape[0], file=f)
        print("There are a total of %d eQTL shared between ancestry" %
              len(shared), file=f)


if __name__ == "__main__":
    main()
