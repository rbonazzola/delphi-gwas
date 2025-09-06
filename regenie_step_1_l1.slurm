#!/bin/bash
#SBATCH -J regenie_step_1_l1 
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=0-119
#SBATCH -o logs/step1_l1/%A_%a.out

set -euo pipefail

EMBEDDING_SIZE=120
BSIZE=${BSIZE:=500}

BFILE=/nfs/research/birney/projects/association/snp_gwas/regenie/resources/ukb22828_allChr_b0_v3_maf01_04_merge
PHENOFILE=pheno_excluding_rel/embeddings_${AGE}_excl_rel_plink.tsv
COVARIATES_FILE=data/covariates_for_gwas_plink.tsv

PHENOTYPE=$(printf "embedding_%03d\n" $SLURM_ARRAY_TASK_ID)

OUTDIR="split_jobs_stage1_l0/age${AGE}/${PHENOTYPE}"

PCS=$(seq 1 20 | sed 's/^/pc/' | tr "\n" ",")
COVARIATES=${PCS},sex

OFILEPAT="loco/emb${EMBEDDING_SIZE}_${AGE}_${PHENOTYPE}"
if [ -f "${OFILEPAT}_1.loco" ]; then
    echo "Error: File '$FILE' already exists. Exiting."
    exit 1
fi

regenie \
  --bed $BFILE \
  --covarFile $COVARIATES_FILE \
  --covarColList $COVARIATES \
  --phenoFile $PHENOFILE \
  --phenoCol $PHENOTYPE \
  --bsize ${BSIZE} \
  --step 1 \
  --keep-l0 \
  --apply-rint \
  --out ${OFILEPAT} \
  --lowmem \
  --lowmem-prefix tmp_rg \
  --run-l1 ${OUTDIR}/split.master
