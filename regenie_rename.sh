#!/bin/bash
AGE=20
EMBEDDING_SIZE=120

BASE="split_jobs_stage1_l0"
OUTDIR="${BASE}/age${AGE}"
mkdir -p "$OUTDIR"

for HALF in head tail; do
  INDIR="split_jobs_age${AGE}_emb${EMBEDDING_SIZE}_stage1_l0_${HALF}"
  for f in ${INDIR}/split_job*_l0_Y*; do
    # extraer el número del fenotipo al final del nombre YN
    Y=$(echo "$f" | sed -E 's/.*_Y([0-9]+)$/\1/')
    if [ "$HALF" == "head" ]; then
      PHENO=$((Y - 1))        # head: PHENOTYPE - 1
    else
      PHENO=$((Y + 59))       # tail: PHENOTYPE + 59
    fi

    # formatear con ceros a la izquierda (ej. 000, 001, ...)
    PHENO_NAME=$(printf "embedding_%03d" "$PHENO")

    # crear subdirectorio por fenotipo
    mkdir -p "${OUTDIR}/${PHENO_NAME}"

    # nombre destino sin el sufijo _Y...
    NEWF=$(echo "$f" | sed -E 's|.*/||; s/_Y[0-9]+$//')

    # crear symlink en vez de mover
    CMD="ln -s $PWD/$f ${OUTDIR}/${PHENO_NAME}/${NEWF}_Y1"
    echo $CMD
    $CMD
    # echo "ln -s $PWD/$f -> ${OUTDIR}/${PHENO_NAME}/${NEWF}_Y1"
  done
done

