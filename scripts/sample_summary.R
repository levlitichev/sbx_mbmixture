suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tools))

# check if script is being run via snakemake
if (!exists("snakemake"))
    stop(sprintf("This script only works as part of a Snakemake pipeline.")) 

readcounts_path <- snakemake@input[["readcounts_path"]]
imp_snp_path <- snakemake@input[["imp_snp_path"]]
sample <- snakemake@params[["sample"]]
chr <- as.integer(snakemake@params[["chr"]])
out_path <- snakemake@output[[1]]
threads <- as.integer(snakemake@threads)

# make sure read counts file exists
if (!file.exists(readcounts_path))
    stop(sprintf("Path to read counts doesn't exist. readcounts_path: %s", readcounts_path))

# make sure imputed SNP file exists
if (!file.exists(imp_snp_path))
    stop(sprintf("Path to imputed SNPs doesn't exist. imp_snp_path: %s", imp_snp_path))

# make sure chromosome was provided
if (is.null(chr))
    stop("Must provide a chromosome.")
cat("Chromosome", chr, "\n")

# make sure output path can be written to
if (!dir.exists(dirname(out_path)))
    stop(sprintf("Can't write to the output path. out_path: %s", out_path))

# import read counts
counts.df <- read.csv(readcounts_path,
                      colClasses=list(chr="character")) # keep the chromosome column as character type
cat("Loaded", readcounts_path, "\n")

# import imp_snps (slow)
imp_snps <- read_csv(imp_snp_path, col_types=cols(chr="i", pos="i", .default="d"), num_threads=threads)
mouse_IDs <- colnames(imp_snps)[3:ncol(imp_snps)]
cat("Loaded", imp_snp_path, "\n")

# SNPs in counts.df should be a subset of imputed bi-allelic SNPs
# (because in get_readcounts.R, we already subset to bi-allelic SNPs in the CC sqlite db)
merged  <- merge(imp_snps, counts.df, by=c("chr", "pos"))
stopifnot(nrow(merged) > 0)
cat("Using", nrow(merged), "SNPs to compare this microbiome sample to genotyped mice...\n")

# sample_results compares this microbiome sample to each genotyped mouse
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

# initialize output
sample_results <- array(0, dim=c(length(mouse_IDs), 3, 2))
dimnames(sample_results) <- list(mouse_IDs, c("AA", "AB", "BB"), c("A", "B"))

# loop over mice
for (ii in 1:length(mouse_IDs)) {
  cat("Comparing microbiome sample to mouse", ii, "of", length(mouse_IDs), "\n") 
  this.mouse <- mouse_IDs[ii]
  
  num.A.for.AA <- sum(merged[merged[[this.mouse]] == 0, "count1"])
  num.B.for.AA <- sum(merged[merged[[this.mouse]] == 0, "count2"])
  num.A.for.AB <- sum(merged[merged[[this.mouse]] == 0.5, "count1"])
  num.B.for.AB <- sum(merged[merged[[this.mouse]] == 0.5, "count2"])
  num.A.for.BB <- sum(merged[merged[[this.mouse]] == 1, "count1"])
  num.B.for.BB <- sum(merged[merged[[this.mouse]] == 1, "count2"])
  sample.result.for.ii <- matrix(c(
    num.A.for.AA, num.B.for.AA,
    num.A.for.AB, num.B.for.AB,
    num.A.for.BB, num.B.for.BB),
    ncol=2, byrow=T)
  sample_results[ii,,] <- sample.result.for.ii
}

# write output
saveRDS(sample_results, out_path)
cat("Saved sample_results for chr", chr, "\n")

