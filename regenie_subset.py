#!/usr/bin/env python3
import sys
import os
import pandas as pd
from tqdm import tqdm

bad_lines = []
def log_bad_line(line):
    print("Bad line:", line)   # la imprime
    bad_lines.append(line)     # la guarda en una lista
    return None 

if len(sys.argv) != 3:
    sys.stderr.write("Usage: filter_regenie_snps.py <embedding_dim> <age>\n")
    sys.exit(1)

embedding_dim = int(sys.argv[1])
age = int(sys.argv[2])

indir = os.path.expandvars("$NB/merged_retry")
outdir = os.path.expandvars("$NB/merged_retry/subsetted")
snplist = os.path.expandvars("$HOME/Delphi/gwas/data/signif_snps.txt")
snplist = os.path.expandvars("$HOME/signif_snps2.txt")

os.makedirs(outdir, exist_ok=True)

# File names
pheno = f"embedding_{embedding_dim:03d}_{age}"
infile = os.path.join(indir, f"{pheno}.regenie")
outfile = os.path.join(outdir, f"{pheno}.filtered")

# Read list of SNPs
snps = set(pd.read_csv(snplist, header=None, sep="\s+").iloc[:,2])

chunksize = 10**6

dfs = []
# for chunk in tqdm(pd.read_csv(infile, sep="\s+", chunksize=chunksize, on_bad_lines=log_bad_line, engine="python"), desc="Loading"):
for chunk in tqdm(pd.read_csv(infile, sep="\s+", chunksize=chunksize, on_bad_lines='warn', engine="c"), desc="Loading"):
    dfs.append(chunk)
df = pd.concat(dfs, ignore_index=True)

# Read REGENIE output
# df = pd.read_csv(infile, sep="\s+", on_bad_lines=log_bad_line, engine="python")
# print(f"Reading: {infile=}")
# import ipdb; ipdb.set_trace()

# Filter
df_filt = df[df["ID"].isin(snps)].copy().assign(embedding=f"embedding_{embedding_dim:03d}", age=age)

# Save
df_filt.to_csv(outfile, sep="\t", index=False)
print(f"Filtered {len(df_filt)} SNPs for {pheno} → {outfile}")
