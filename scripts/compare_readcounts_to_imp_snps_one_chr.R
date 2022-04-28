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
sample_results_out_path <- snakemake@output[["sample_results_out_path"]]
pair_results_out_path <- snakemake@output[["pair_results_out_path"]]

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

# make sure output paths can be written to
if (!dir.exists(dirname(sample_results_out_path)))
    stop(sprintf("Can't write to the sample_results output path. sample_results_out_path: %s", sample_results_out_path))
if (!dir.exists(dirname(pair_results_out_path)))
    stop(sprintf("Can't write to the pair_results output path. pair_results_out_path: %s", pair_results_out_path))

# import read counts
counts.df <- read.csv(readcounts_path,
                      colClasses=list(chr="character")) # keep the chromosome column as character type
cat("Loaded", readcounts_path, "\n")

# the sample should be named like `DO_1D_3011_044w`
# and the corresponding mouse ID will be `DO-1D-3011`
mouse.ID.for.this.mb.sample <- paste0(strsplit(sample, "_")[[1]][1:3], collapse="-")
cat("Mouse ID for this microbiome sample:", mouse.ID.for.this.mb.sample, "\n")

# import imp_snps (slow)
imp_snps <- read_csv(imp_snp_path, col_types=cols(chr="i", pos="i", .default="d"))
mouse_IDs <- colnames(imp_snps)[3:ncol(imp_snps)]
cat("Loaded", imp_snp_path, "\n")

# check that the mouse ID for our microbiome sample matches a mouse ID in imp_snps
if (!(mouse.ID.for.this.mb.sample %in% mouse_IDs))
  stop(sprintf("No mouse ID in imp_snps matches our mouse ID. Mouse ID in microbiome sample: %s, first 3 mice in imp_snps: %s", mouse.ID.for.this.mb.sample, paste0(mouse_IDs[1:3], collapse=", ")))

# SNPs in counts.df should be a subset of imputed bi-allelic SNPs
# (because in get_readcounts.R, we already subset to bi-allelic SNPs in the CC sqlite db)
merged  <- merge(imp_snps, counts.df, by=c("chr", "pos"))
stopifnot(nrow(merged) > 0)
cat("Using", nrow(merged), "SNPs to compare this microbiome sample to genotyped mice...\n")

# initialize outputs
# sample_results compares this microbiome sample to each genotyped mouse
sample_results <- array(0, dim=c(length(mouse_IDs), 3, 2))
dimnames(sample_results) <- list(mouse_IDs, c("AA", "AB", "BB"), c("A", "B"))

# pair_results is used for determining whether one sample is a mixture of samples
# for all genotypes, simultaneously gets counts for the expected mouse genotype versus another genotype
pair_results <- array(0, dim=c(length(mouse_IDs), 3, 3, 2))
dimnames(pair_results) <- list(mouse_IDs, c("AA", "AB", "BB"),
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
for (ii in 1:length(mouse_IDs)) {
  cat("Comparing microbiome sample to mouse", ii, "of", length(mouse_IDs), "\n") 
  this.mouse <- mouse_IDs[ii]
  
  # SAMPLE_RESULTS
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
  
  # PAIR_RESULTS
  # count up the number of reads for major and minor alleles for every combination of
  # expected genotype versus all other genotypes
  pair.result.for.ii <- merged %>%
    mutate(this.genotype=factor(.data[[this.mouse]], levels=c(0,0.5,1)),
           expected.genotype=factor(.data[[mouse.ID.for.this.mb.sample]], levels=c(0,0.5,1))) %>%
    dplyr::filter(!is.na(this.genotype)) %>%
    dplyr::filter(!is.na(expected.genotype)) %>%
    group_by(this.genotype, expected.genotype, .drop=F) %>%
    summarise(A=sum(count1), B=sum(count2), .groups="drop") 
  pair_results[ii,,,1] <- pair.result.for.ii  %>%
    pivot_wider(id_cols=expected.genotype, names_from=this.genotype, values_from=A) %>% 
    dplyr::select(c(`0`,`0.5`,`1`)) %>% mutate_all(function(x) as.integer(x)) %>% as.matrix
  pair_results[ii,,,2] <- pair.result.for.ii  %>%
    pivot_wider(id_cols=expected.genotype, names_from=this.genotype, values_from=B) %>%
    dplyr::select(c(`0`,`0.5`,`1`)) %>% mutate_all(function(x) as.integer(x)) %>% as.matrix
}  

# write output for this chromosome  
saveRDS(sample_results, sample_results_out_path)
saveRDS(pair_results, pair_results_out_path)
cat("Saved results for chr", chr, "\n")
