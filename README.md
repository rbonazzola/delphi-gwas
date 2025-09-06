# Delphi GWAS
Scripts to perform GWAS on Delphi internal representation using `regenie`.

# Steps

## Preliminaries
  - Generate phenotype files: `PHENOFILES=...`
  - Generate covariates file: `COVARIATES_FILE=...` (must include `FID` and `IID` columns)
  - Generate include and exclude list of individuals: `INCLUDE_IDS=...`, `EXCLUDE_IDS=...`
  - Locate relatedness file: `RELFILE=...`
  - Locate `BFILE=...`
  - Locate `BGEN=...` and `BGEN_SAMPLE=...`
  - Install regenie and put in `PATH`.

## Preprocessing:
  - Define kinship threshold.
  - Run subject filtering step.
  - Compile phenotypes into a single file. This may require adding suffixes in case of name collision (optional)

## `regenie: **Step 1 / level 0**:
  - Format phenotype file: add `FID` and `IID`.
  - Decide `BSIZE=...` for step 1/level 0.
  - Define rule for phenotype names.
  - Define output directory.


## `regenie`: **Step 1 / level 1**
  - "Linkify" files: put them into a different folder structure, where folder contains the phenotype index, and the original phenotype suffix is changed into `_Y1`.
  - Modify master to reflect this change.

## `regenie`: **Step 2** (association tests)
  - Generate regions file. Define region length. Filter out regions with no variants.
  - Define `minMAC` and `minINFO`.
  - Choose `NTHREADS=...`
  - Define `BSIZE` (determines memory usage)
  - Execute step 2 (as job array)
  - Detect failed jobs and re-run.

## Post-processing
  - Merge results for different regions (one file per phenotype)
  - Find significant SNPs (SNPs that are genome-wide significant for at least one embedding dimension and age)
  - Filter results for the previous SNPs and compile them into a single file, one file per (SNP, age) and one column per embedding dimension (R script).
