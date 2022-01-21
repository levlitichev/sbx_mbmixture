suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tools))

# check if script is being run via snakemake
if (!exists("snakemake"))
    stop(sprintf("This script only works as part of a Snakemake pipeline.")) 

args <- list()
args$readcounts_path <- snakemake@input[["readcounts_path"]]
args$chr <- chr <- as.integer(snakemake@params[["chr"]])
args$sample_results_out_path <- snakemake@output[["sample_results_out_path"]]
args$pair_results_out_path <- snakemake@output[["pair_results_out_path"]]
args$snp_db_path <- snakemake@config[["sbx_mbmixture"]][["CC_sqlite_db"]]
args$imp_snp_dir <- snakemake@config[["sbx_mbmixture"]][["imp_snps_dir"]]

# make sure read counts file exists
if (!file.exists(args$readcounts_path))
    stop(sprintf("Path to read counts doesn't exist. args$readcounts_path: %s", args$readcounts_path))

# make sure imputed SNPs directory exists
if (!dir.exists(args$imp_snp_dir))
    stop(sprintf("Can't find directory with imputed SNPs. args$imp_snp_dir: %s", args$imp_snp_dir))

# make sure CC SNP db exists
if (!file.exists(args$snp_db_path))
    stop(sprintf("Path to CC SNP db doesn't exist. args$snp_db_path: %s", args$snp_db_path))

# make sure output paths can be written to
if (!dir.exists(dirname(args$sample_results_out_path)))
    stop(sprintf("Can't write to the sample_results output path. args$sample_results_out_path: %s", args$sample_results_out_path))
if (!dir.exists(dirname(args$pair_results_out_path)))
    stop(sprintf("Can't write to the pair_results output path. args$pair_results_out_path: %s", args$pair_results_out_path))

# the readcounts file should be named like `DO_1D_3011_044w.csv`
# and the corresponding mouse ID is `DO-1D-3011`
mouse.ID.for.this.mb.sample <- paste0(strsplit(basename(args$readcounts_path), "_")[[1]][1:3], collapse="-")
cat("Mouse ID for this microbiome sample:", mouse.ID.for.this.mb.sample, "\n")

# import read counts
counts.df.all <- read.csv(args$readcounts_path,
                 colClasses=list(chr="character")) # keep the chromosome column as character type
cat("Loaded", args$readcounts_path, "\n")

# connect to CC database
snp.db <- dbConnect(SQLite(), args$snp_db_path)

# subset read counts to just this chromosome
counts.df <- counts.df.all[counts.df.all$chr == chr, ]  

# extract positions and IDs for CC SNPs 
cat("Loading CC SNPs...")
CC_snp_info <- dbGetQuery(snp.db, paste0("select pos,snp_id from variants where chr=='", chr, "'"))
cat(" done.\n")

# add snp_id to counts.df
counts.df.w.SNP.id <- merge(counts.df, CC_snp_info, by="pos")

# some positions correspond to more than 1 snp_id; remove duplicate snp_ids
counts.df.w.SNP.id <- counts.df.w.SNP.id[!duplicated(counts.df.w.SNP.id$pos), ]
stopifnot(length(unique(counts.df.w.SNP.id$pos)) == length(counts.df.w.SNP.id$pos))
stopifnot(nrow(counts.df.w.SNP.id) == nrow(counts.df))

# make sure imputed SNPs file exists for this chromosome
imp_snp_path <- file.path(args$imp_snp_dir, paste0("imp_snp_", chr, ".csv"))
if (!file.exists(imp_snp_path))
  stop(sprintf("Can't find imputed SNPs for chromosome %s. imp_snp_path: %s", chr, imp_snp_path))

# import imp_snps
cat("Loading imputed SNPs...")
imp_snps <- read.csv(file.path(args$imp_snp_dir, paste0("imp_snp_", chr, ".csv")), row.names=1, check.names=F)
cat(" done.\n")

# check that the mouse ID for our microbiome sample matches a mouse ID in imp_snps
if (!(mouse.ID.for.this.mb.sample %in% rownames(imp_snps)))
  stop(sprintf("No mouse ID in imp_snps matches our mouse ID. First 3 mice in imp_snps: %s", paste0(rownames(imp_snps)[1:3], collapse=", ")))

# subset to SNPs that we were able to impute, i.e. in colnames(imp_snps)
SNP.ids <- intersect(counts.df.w.SNP.id$snp_id, colnames(imp_snps))
stopifnot(length(SNP.ids) > 0)
counts.df.w.SNP.id <- counts.df.w.SNP.id[match(SNP.ids, counts.df.w.SNP.id$snp_id), ]
imp_snps <- imp_snps[, SNP.ids]
stopifnot(all(counts.df.w.SNP.id$snp_id == colnames(imp_snps)))
cat("Using", length(SNP.ids), "SNPs to compare this microbiome sample to genotyped mice...\n")

# initialize outputs
# sample_results compares this microbiome sample to each genotyped mouse
sample_results <- array(0, dim=c(nrow(imp_snps), 3, 2))
dimnames(sample_results) <- list(rownames(imp_snps), c("AA", "AB", "BB"), c("A", "B"))

# pair_results is used for determining whether one sample is a mixture of samples
# for all genotypes, simultaneously gets counts for the expected mouse genotype versus another genotype
pair_results <- array(0, dim=c(nrow(imp_snps), 3, 3, 2))
dimnames(pair_results) <- list(rownames(imp_snps), c("AA", "AB", "BB"),
                               c("AA", "AB", "BB"), c("A", "B"))

# example of how sample_results will look
# sample_results[1,,]:
#         A     B
# AA   5877   646
# AB   2835  1145
# BB    398   340
# how to interpret: 
# for the 1st mouse, we found:
#     5877 reads with the major allele (A AKA count1) at SNPs that were homozygous for the major allele (AA) (these are 'correct' reads);
#     646 reads with the minor allele (B AKA count2) across these same AA SNPs ('incorrect' reads);
#     2835 reads with the major allele A at heterozygous SNPs
#     (we don't really know what to do with this, so we will ignore these counts);
#     398 reads with the major allele A at SNPs homozygous for the minor allele (BB) ('incorrect' reads); and
#     340 reads with the minor allele B at these same BB SNPs ('correct' reads)

# loop over mice
for (ii in 1:nrow(imp_snps)) {
  
  # SAMPLE_RESULTS
  # add this mouse's genotype to counts.df
  this.genotype.df <- t(imp_snps[ii, ])
  colnames(this.genotype.df) <- "this.genotype"
  counts.df.for.ii <- merge(
    counts.df.w.SNP.id, this.genotype.df,
    by.x="snp_id", by.y="row.names")
  
  # count up the number of reads for major (A) and minor (B) alleles
  # separately for homozygous major (AA), heterozygous (AB), and homozygous minor (BB) SNPs
  sample.result.for.ii <- counts.df.for.ii %>% 
    group_by(this.genotype) %>%
    summarise(A=sum(count1), B=sum(count2), .groups="drop") %>% as.matrix
  if (any(dim(sample.result.for.ii) != c(3,2)))
    stop(sprintf("Didn't find all genotypes (AA, AB, BB). Probably insufficient reads for this sample."))
  sample_results[ii,,] <- sample.result.for.ii[, c("A","B")]
  
  # PAIR_RESULTS
  # add the expected mouse genotype
  expected.mouse.genotype.df <- t(imp_snps[mouse.ID.for.this.mb.sample, ])
  colnames(expected.mouse.genotype.df) <- "expected.genotype"
  counts.df.for.ii <- merge(
    counts.df.for.ii, expected.mouse.genotype.df,
    by.x="snp_id", by.y="row.names")
  
  # count up the number of reads for major and minor alleles for every combination of
  # expected genotype versus all other genotypes
  pair.result.for.ii <- counts.df.for.ii %>%
    group_by(this.genotype, expected.genotype) %>%
    summarise(A=sum(count1), B=sum(count2), .groups="drop") 
  pair_results[ii,,,1] <- pair.result.for.ii  %>%
    pivot_wider(id_cols=expected.genotype, names_from=this.genotype, values_from=A) %>% 
    dplyr::select(c(`1`,`2`,`3`)) %>% as.matrix
  pair_results[ii,,,2] <- pair.result.for.ii  %>%
    pivot_wider(id_cols=expected.genotype, names_from=this.genotype, values_from=B) %>%
    dplyr::select(c(`1`,`2`,`3`)) %>% as.matrix
}  

# write output for this chromosome  
saveRDS(sample_results, args$sample_results_out_path)
saveRDS(pair_results, args$pair_results_out_path)
cat("Saved results for chr", chr, ".\n")

# disconnect from CC database
dbDisconnect(snp.db)
