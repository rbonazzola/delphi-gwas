#!/bin/bash
AGES="20 30 40 50 60"          
EMBEDDING_SIZE=120

BASE="split_jobs_stage1_l0"

for AGE in $AGES; do
  for PHENO in $(seq 0 119); do

      if [ $PHENO -lt 60 ]; then HALF=head; else HALF=tail; fi

      PHENONAME=$(printf "embedding_%03d" "$PHENO")

      IMASTER="${PWD}/split_jobs_age${AGE}_emb${EMBEDDING_SIZE}_stage1_l0_${HALF}/split.master"
      SRCDIR=$(dirname $IMASTER)
      SRCDIR=$(basename $SRCDIR)

      OUTDIR="${BASE}/age${AGE}/${PHENONAME}"

      OMASTER="${OUTDIR}/split.master"
      DSTDIR=$(dirname $OMASTER)
      DSTDIR=$(echo $DSTDIR | cut -d'/' -f 1,2,3)

      cp "${IMASTER}" "${OMASTER}"
      sed -i "s|${SRCDIR}|${DSTDIR}|g" ${OMASTER} 
  done
done
