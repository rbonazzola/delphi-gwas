#!/usr/bin/env bash
#SBATCH -J merge_pheno
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH --array=1-600
#SBATCH -o logs/merged_retry/%a.out

NOBACKUP=$NB
cd $NOBACKUP

REGIONS_FILE=${HOME}/Delphi/gwas/data/regions_2mb_hg19.bed
OUTDIR=${NOBACKUP}/merged_retry
mkdir -p $OUTDIR

# Generate phenotype list from header (columns 3–602)
PHENOS=$(head -n1 /homes/bonazzola/Delphi/gwas/pheno_excluding_rel/merged.tsv | cut -f3-602 | tr '\t' '\n')
PHENO=$(echo "$PHENOS" | sed -n "${SLURM_ARRAY_TASK_ID}p")

OUTFILE=${OUTDIR}/${PHENO}.regenie
FAILFILE=${OUTDIR}/${PHENO}.failed_regions.txt
TMPFILE=$(mktemp)

IDIR=${NOBACKUP}/emb120_regenie2

# Header from the first region file
read chr start end < <(head -n1 $REGIONS_FILE)
FIRSTFILE=${IDIR}/emb120_chr${chr}_${start}-${end}_${PHENO}.regenie
head -n1 "$FIRSTFILE" > $TMPFILE

# Reset failed region list
> $FAILFILE

# Concatenate all regions
while read chr start end; do
    FILE=${IDIR}/emb120_chr${chr}_${start}-${end}_${PHENO}.regenie
    if [[ -f $FILE ]]; then
        tail -n +2 $FILE >> $TMPFILE
    else
        echo -e "${chr}\t${start}\t${end}" >> $FAILFILE
    fi
done < $REGIONS_FILE

mv $TMPFILE $OUTFILE

# Report outcome
if [[ ! -s $FAILFILE ]]; then
    rm -f $FAILFILE
    echo "Merged $PHENO → $OUTFILE (no missing regions)"
else
    echo "Merged $PHENO → $OUTFILE (missing $(wc -l < $FAILFILE) regions, see $FAILFILE)"
fi

awk '$13 > 7.3' $OUTFILE > ${OUTFILE%*.regenie}_signif.regenie
