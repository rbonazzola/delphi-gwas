#!/bin/bash
#SBATCH -J regenie_stage_2
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --array=1-5
#SBATCH -o logs/step2_retry/age20_60/%A_%a.out

set -euo pipefail

EMBEDDING_SIZE=120

i=$SLURM_ARRAY_TASK_ID
read CHROMOSOME START END < <(sed -n "${i}p" data/regions_2mb_hg19.bed)
REGION=${CHROMOSOME}:${START}-${END}

BFILE=/nfs/research/birney/projects/association/snp_gwas/regenie/resources/ukb22828_allChr_b0_v3_maf01_04_merge
BGEN=/nfs/research/sds/sds-ukb-geno/genotyping/imputed/ukb22828_c${CHROMOSOME}_b0_v3.bgen
NOBACKUP=/hps/nobackup/birney/users/bonazzola

echo $REGION

COVARIATES_FILE=data/covariates_for_gwas_plink.tsv

OFILEPAT="${NOBACKUP}/emb120_regenie2/emb${EMBEDDING_SIZE}_chr${CHROMOSOME}_${START}-${END}"
PHENOFILE="${NOBACKUP}/pheno_excluding_rel/merged.tsv"
PHENOTYPES=$(head -n1 ${PHENOFILE} | cut -f3-602 | sed 's/\t/,/g')

echo $OFILEPAT

CMD="regenie \
  --bgen $BGEN \
  --sample /nfs/research/sds/sds-ukb-geno/genotyping/imputed/sample/ukb22828_c1_b0_v3_s487276.sample \
  --covarFile $COVARIATES_FILE \
  --covarColList pc{1:20},sex \
  --phenoFile $PHENOFILE \
  --phenoColList $PHENOTYPES \
  --range ${REGION} \
  --step 2 \
  --apply-rint \
  --out ${OFILEPAT} \
  --minINFO 0.3 \
  --minMAC 6000 \
  --bsize 2000 \
  --threads 48 \
  --pred ${NOBACKUP}/loco/all_pred_with_age.list"

$CMD
