#!/bin/bash
#SBATCH --job-name=linkify
#SBATCH --array=0-119           # 120 fenotipos por edad
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -o logs/linkify/%A_%a.out
#SBATCH -e logs/linkify/%A_%a.err

AGE=$1
EMBEDDING_SIZE=120

BASE="split_jobs_stage1_l0"
OUTDIR="${BASE}/age${AGE}"
mkdir -p "$OUTDIR"

# distinguir head/tail según el índice
if [ $SLURM_ARRAY_TASK_ID -lt 60 ]; then
    HALF=head
    Y=$((SLURM_ARRAY_TASK_ID+1))   # regenie arranca en Y1
    PHENO=$((Y-1))                 # 0–59
else
    HALF=tail
    Y=$((SLURM_ARRAY_TASK_ID-59))  # también arranca en Y1
    PHENO=$((Y+59))                # 60–119
fi

PHENO_NAME=$(printf "embedding_%03d" "$PHENO")

INDIR="split_jobs_age${AGE}_emb${EMBEDDING_SIZE}_stage1_l0_${HALF}"
mkdir -p "${OUTDIR}/${PHENO_NAME}"

for f in ${INDIR}/split_job*_l0_Y${Y}; do
    NEWF=$(basename "$f" | sed -E 's/_Y[0-9]+$//')
    # rm $NEWF
    ln -s "$PWD/$f" "${OUTDIR}/${PHENO_NAME}/${NEWF}_Y1"
done
