import pandas as pd
import sys
phenotype_file = sys.argv[1]

df = pd.read_csv(phenotype_file, sep='\t')
if len(df.columns) == 1:
    df = pd.read_csv(phenotype_file)

df.insert(0, "FID", df["ID"])
df.rename(columns={"ID": "IID"}, inplace=True)

df.to_csv(phenotype_file.replace(".csv", "")+"_plink.tsv", sep="\t", index=False, na_rep="NA")
