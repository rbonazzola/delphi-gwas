#!/bin/bash
#SBATCH -J widen_gwas
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH -o logs/widen_%A_%a.out
#SBATCH -e logs/widen_%A_%a.err

set -euo pipefail

BFILE=/nfs/research/birney/projects/association/snp_gwas/regenie/resources/ukb22828_allChr_b0_v3_maf01_04_merge
PHENO=pheno_excluding_rel/embeddings_20_excl_rel.csv
PHENO=pheno_excluding_rel/embeddings_20_excl_rel_plink.tsv
COVARIATES=data/covariates_for_gwas.csv
COVARIATES=data/covariates_for_gwas_plink.tsv

regenie \
  --step 1 \
  --bed $BFILE \
  --covarFile $COVARIATES \
  --phenoFile $PHENO \
  --phenoCol embedding_000 \
  --bsize 100 \
  --out fit_bin_000 \
  --lowmem \
  --lowmem-prefix tmp_rg_000
  #--remove example/fid_iid_to_remove.txt \
