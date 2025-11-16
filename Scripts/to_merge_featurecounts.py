#conda install pandas
import pandas as pd
import glob

files = glob.glob("*.tabular")

dfs = []
for f in files:
    df = pd.read_csv(f, sep="\t")
    dfs.append(df)

merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on="Geneid")

merged.to_csv("gene_counts_matrix.csv", index=False)

