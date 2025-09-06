AGE=$1

python remove_related.py \
  --phenotype_file /home/bonazzola/Delphi/embeddings/data/embedding_120/embeddings_${AGE}.csv \
  --relatedness_file /home/bonazzola/Delphi/data/geno/ukb_rel_a49978_s488237.dat \
  --bgen_sample_file /home/bonazzola/Delphi/data/geno/sample/ukb22828_c1_b0_v3_s487276.sample \
  --keep_file /home/bonazzola/Delphi/data/datasets/genetic_white_ids.txt \
  --na_code NA \
  --output_dir /home/bonazzola/Delphi/gwas/pheno_excluding_rel \
  --tmpdir /home/bonazzola/Delphi/gwas/tmp_pheno \
  --seed 42 \
  --gzip_output \
  --output_file_prefix embeddings_${AGE}_excl_rel
