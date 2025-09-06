#!/bin/bash
#SBATCH -J regenie_stage_1 
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1-1000
#SBATCH -o logs/step1_l0/120_20_%A_%a.out
#SBATCH -e logs/step1_l0/120_20_%A_%a.err

set -euo pipefail

EMBEDDING_SIZE=120
AGE=${AGE:=20}
BSIZE=${BSIZE:=500}

BFILE=/nfs/research/birney/projects/association/snp_gwas/regenie/resources/ukb22828_allChr_b0_v3_maf01_04_merge
PHENOFILE=pheno_excluding_rel/embeddings_${AGE}_excl_rel_plink.tsv
COVARIATES_FILE=data/covariates_for_gwas_plink.tsv

WHICH_HALF=${WHICH_HALF:=head}
PHENOTYPES=$(seq -w 0 119 | $WHICH_HALF -n60 | sed 's/^/embedding_/' | tr '\n' ',')

if [ -n "${SLURM_ARRAY_TASK_ID-}" ]; then
  echo "ANDAAA"
fi

echo $WHICH_HALF
echo $AGE

OUTDIR="split_jobs_age${AGE}_emb${EMBEDDING_SIZE}_stage1_l0_${WHICH_HALF}"
PREFIX=$OUTDIR/split

mkdir -p $OUTDIR

PCS=$(seq 1 20 | sed 's/^/pc/' | tr "\n" ",")
COVARIATES=${PCS},sex

if [ -n "${SLURM_ARRAY_TASK_ID-}" ]; then
  JOBSPEC_OR_EXEC="--run-l0 ${PREFIX}.master,${SLURM_ARRAY_TASK_ID}"
else
  JOBSPEC_OR_EXEC="--split-l0 ${PREFIX},1000"
fi

regenie \
  --bed $BFILE \
  --covarFile $COVARIATES_FILE \
  --covarColList $COVARIATES \
  --phenoFile $PHENOFILE \
  --phenoColList $PHENOTYPES \
  --bsize ${BSIZE} \
  --step 1 \
  --apply-rint \
  --out emb${EMBEDDING_SIZE}_${AGE}_test_split \
  --lowmem \
  --lowmem-prefix tmp_rg \
  $JOBSPEC_OR_EXEC
