# %%
import os, sys
os.getcwd()
sys.path.append("gwas")
import pandas as pd

from importlib import reload 
import gwas_covariates_helpers

reload(gwas_covariates_helpers)
gcov = gwas_covariates_helpers

# %%
import yaml

yaml_str = f'''
../data/datasets/genetic_pcs_22009.txt:
  - id: "f.eid"
''' \
+ "\n".join([f'  - f.22009.0.{i}: pc{i}' for i in range(1, 41)]) + \
'''
../data/datasets/sex_31.txt:
  - id: "f.eid"
  - f.31.0.0: sex

../data/datasets/height.csv:
  - id : "f.eid"
  - reduce:
      columns: [f.50.0.0, f.50.1.0, f.50.2.0, f.50.3.0]
      name: height
      method: mean
'''

# %%
print(yaml_str)

# %%
covariates = yaml.safe_load(yaml_str)
print(covariates)

cov_df = gcov.generate_covariates_df(covariates_config=covariates, return_individual_dfs=False)
# cov_df, all_dfs = gcov.generate_covariates_df(covariates_config=covariates, return_individual_dfs=True)

# %%
os.getcwd()
# %%
ethn_df = pd.read_csv("~/Delphi/data/datasets/22006.csv")
ethn_df = ethn_df.rename({"eid": "ID", "22006-0.0": "white"}, axis=1)
ethn_df.white = ethn_df.white.notnull().astype(int)

emb120_df = pd.read_csv("/home/bonazzola/Delphi/embeddings/data/embedding_120/embeddings_20.csv")
emb120_df

# %%
all_df = pd.concat([
    emb120_df.rename({"subject_id": "ID"}, axis=1).set_index("ID").rename_axis(None).rename(index=lambda x: int(x)),
    cov_df.set_index("ID").rename_axis(None).rename(index=lambda x: int(x)),
    ethn_df.set_index("ID").rename_axis(None).rename(index=lambda x: int(x))
], axis=1)

( all_df := all_df.query("white == 1").drop("white", axis=1) ).sample(10)
# %%
covariates = (all_pcs := [f'pc{i}' for i in range(1, 41)]) + ['sex']
covariates

covariates_df = all_df[covariates]
covariates_df.to_csv("covariates_for_gwas.csv", index=True)

# %%
from tqdm import tqdm

covariates = (all_pcs := [f'pc{i}' for i in range(1, 41)]) + ['sex']
embedding_col = lambda i: f'embedding_{str(i).zfill(3)}'
import gc

fit = {}
fit_summaries = []
residues = []

for j in tqdm(range(0, 120)):
    model = gcov.fit_linear_model(all_df, embedding_col(j), covariates)
    fit_results_df = pd.DataFrame({
      "coef": model.params,
      "std_err": model.bse,
      "t": model.tvalues,
      "pval": model.pvalues,
    })
    fit_summaries.append(fit_results_df)
    residues.append( model.resid )

    gc.collect()

# %%
fit_summary_df = pd.concat([df.assign(embedding=i) for i, df in enumerate(fit_summaries)], axis=0)

# %%
pc_only_pvals = fit_summary_df[['pval', 'embedding']].pivot(columns="embedding", values="pval").drop(['const', 'sex']).loc[[f'pc{i}' for i in range(1, 41)]]

# %%
import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(24, 6))
sns.heatmap(
  -np.log10(pc_only_pvals.values).clip(-40, 0)
);

plt.title("log10(p) of linear regression of embedding dimensions on genetic PCs")
plt.xlabel("Embedding dimension")
plt.ylabel("PCs")

# %%

