"""
This script is used to read excel file.
"""
import pandas as pd

def main():
    fn = "https://static-content.springer.com/esm/art%3A10.1186%2F1471-"+\
        "2164-14-606/MediaObjects/12864_2013_5321_MOESM1_ESM.xlsx"
    df = pd.read_excel(fn, sheet_name="ANOVA results")
    df.to_csv("Korostynski2013_mice.csv", sep=",", index=False)


if __name__ == "__main__":
    main()
