library(dplyr)
library(tidyr)

gwas_df = read.csv("/hps/nobackup/birney/users/bonazzola/merged_retry/subsetted/merged_embeddings.filtered", sep = '\t')

phenotypes = paste0("dim", 1:120)

pvals_wide <- gwas_df %>% select(CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, embedding, age, LOG10P) %>% pivot_wider(names_from=embedding, values_from=LOG10P)
 
pvals_wide <- pvals_wide %>% rename(SNP=ID, CHR=CHROM, BP=GENPOS)

colnames(pvals_wide)[8:127] <- phenotypes
pvals_wide = cbind(pvals_wide[1:7], 10**(-pvals_wide[phenotypes]))

write.csv(x = pvals_wide, file = "gwas_only_significant_snps_delphi120_imputed_pvals.csv", quote = FALSE, row.names = FALSE)

########################################################################
beta_wide <- gwas_df %>% select(CHROM, GENPOS, ID, ALLELE0, ALLELE1, A1FREQ, embedding, age, BETA) %>% pivot_wider(names_from=embedding,values_from = BETA) %>% rename(SNP=ID, CHR=CHROM, BP=GENPOS)
  
colnames(beta_wide)[8:127] <- phenotypes

write.csv(beta_wide, "gwas_only_significant_snps_delphi120_imputed_betas.csv", quote = FALSE, row.names = FALSE)
