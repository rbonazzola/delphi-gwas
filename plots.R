suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(logger)
})

#manhattan <- qqman::manhattan
#qq <- qqman::qq
manhattan <- fastman::fastman
qq <- fastman::fastqq

# —————————————————— Reading arguments ——————————————————
args <- commandArgs(trailingOnly = TRUE)
z <- as.integer(args[1])
pheno <- sprintf("embedding_%03d", z-1)
age <- args[2] 

# —————————————————— Loading GWAS ——————————————————
log_info("Loading GWAS file")
file <- glue::glue('/hps/nobackup/birney/users/bonazzola/merged/{pheno}_{age}.regenie')
gwas <- fread(file, select = c("CHROM", "GENPOS", "ID", "LOG10P"))
gwas$P = 10^(-gwas$LOG10P) 
log_info("GWAS file loaded")
browser()
# —————————————————— Manhattan Plot ——————————————————
png(filename=glue::glue("outputs/emb120/figures/manhattan_{pheno}_{age}.png"),width=3000,height=1200,res=200); 

log_info("Starting Manhattan plot")
manhattan(
  gwas, "CHROM", "GENPOS", "P", "ID", 
  col=c("navy", "skyblue"), 
  suggestiveline = FALSE, 
  genomewideline = -log10(5e-8), 
  cex=0.5,
  main = glue::glue("Embedding {z} at age {age}"),
  ylim = c(0, 30)
)

dev.off()
log_info("Manhattan plot finished")

# —————————————————— QQ Plot ——————————————————
chisq <- qchisq(1 - gwas$P, 1)
lambda_gc <- median(chisq) / qchisq(0.5, 1)

log_info("Starting QQ-plot")
png(filename=glue::glue("outputs/emb120/figures/qqplot_{pheno}_{age}.png"),width=1200,height=1200,res=100);
qq(
  gwas$P,
  main = glue::glue("QQ-plot for dim. {z} at age {age}"),
  #main = bquote("QQ-plot " ~ .(z) ~ " at age " ~ .(age) ~
  #              " (" ~ lambda[GC] == .(round(lambda_gc, 3)) ~ ")"),
  cex = 1
)
dev.off()
log_info("QQ-plot finished")

# —————————————————— The end ——————————————————
