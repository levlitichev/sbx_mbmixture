suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tools))

# check if script is being run via snakemake
if (!exists("snakemake"))
    stop(sprintf("This script only works as part of a Snakemake pipeline.")) 

args <- list()
args$readcounts_path <- "~/DO_1D_3011_148w_readcounts.csv.gz")
args$chr <- chr <- 19
args$sample_results_out_path <- "~/DO_1D_3011_148w_chr19_sample_results.Rds"
args$pair_results_out_path <- "~/DO_1D_3011_148w_chr19_pair_results.Rds"
args$imp_snp_path <- "/project/thaisslab/mbmixture/2021-01_snpcalls/with_mouse_ids/chr19_imputed_w_mouse_ids.csv.gz"

#args <- list()
#args$readcounts_path <- snakemake@input[["readcounts_path"]]
#args$chr <- chr <- as.integer(snakemake@params[["chr"]])
#args$sample_results_out_path <- snakemake@output[["sample_results_out_path"]]
#args$pair_results_out_path <- snakemake@output[["pair_results_out_path"]]
#args$snp_db_path <- snakemake@config[["sbx_mbmixture"]][["CC_sqlite_db"]]
#args$imp_snp_dir <- snakemake@config[["sbx_mbmixture"]][["imp_snps_dir"]]

# make sure read counts file exists
if (!file.exists(args$readcounts_path))
    stop(sprintf("Path to read counts doesn't exist. args$readcounts_path: %s", args$readcounts_path))

# make sure imputed SNP file exists
if (!file.exists(args$imp_snp_path))
    stop(sprintf("Path to imputed SNPs doesn't exist. args$imp_snp_path: %s", args$imp_snp_path))

# make sure chromosome was provided
if (is.null(args$chr))
    stop("Must provide a chromosome.")
cat("Chromosome", args$chr, "\n")

# make sure output paths can be written to
if (!dir.exists(dirname(args$sample_results_out_path)))
    stop(sprintf("Can't write to the sample_results output path. args$sample_results_out_path: %s", args$sample_results_out_path))
if (!dir.exists(dirname(args$pair_results_out_path)))
    stop(sprintf("Can't write to the pair_results output path. args$pair_results_out_path: %s", args$pair_results_out_path))

# the readcounts file should be named like `DO_1D_3011_044w_chr19_readcounts.csv`
# and the corresponding mouse ID is `DO-1D-3011`
mouse.ID.for.this.mb.sample <- paste0(strsplit(basename(args$readcounts_path), "_")[[1]][1:3], collapse="-")
cat("Mouse ID for this microbiome sample:", mouse.ID.for.this.mb.sample, "\n")

# import read counts
counts.df <- read.csv(args$readcounts_path,
                      colClasses=list(chr="character")) # keep the chromosome column as character type
cat("Loaded", args$readcounts_path, "\n")

# import imp_snps
cat("Loading imputed SNPs...")
imp_snps <- read.csv(args$imp_snp_path, check.names=F)
mouse_IDs <- colnames(imp_snps)[3:ncol(imp_snps)]
cat(" done.\n")

# check that the mouse ID for our microbiome sample matches a mouse ID in imp_snps
if (!(mouse.ID.for.this.mb.sample %in% mouse_IDs))
  stop(sprintf("No mouse ID in imp_snps matches our mouse ID. First 3 mice in imp_snps: %s", paste0(mouse_IDs[1:3], collapse=", ")))

# counts.df should be a superset of SNPs that could be imputed
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
  # add the expected mouse genotype
  expected.mouse.genotype.df <- t(imp_snps[mouse.ID.for.this.mb.sample, ])
  colnames(expected.mouse.genotype.df) <- "expected.genotype"
  counts.df.for.ii <- merge(
    counts.df.for.ii, expected.mouse.genotype.df,
    by.x="snp_id", by.y="row.names")
  
  # count up the number of reads for major and minor alleles for every combination of
  # expected genotype versus all other genotypes
  pair.result.for.ii <- counts.df.for.ii %>%
    mutate(this.genotype=factor(this.genotype, levels=c(1,2,3)),
           expected.genotype=factor(expected.genotype, levels=c(1,2,3))) %>% 
    group_by(this.genotype, expected.genotype, .drop=F) %>%
    summarise(A=sum(count1), B=sum(count2), .groups="drop") 
  pair_results[ii,,,1] <- pair.result.for.ii  %>%
    pivot_wider(id_cols=expected.genotype, names_from=this.genotype, values_from=A) %>% 
    dplyr::select(c(`1`,`2`,`3`)) %>% mutate_all(function(x) as.integer(x)) %>% as.matrix
  pair_results[ii,,,2] <- pair.result.for.ii  %>%
    pivot_wider(id_cols=expected.genotype, names_from=this.genotype, values_from=B) %>%
    dplyr::select(c(`1`,`2`,`3`)) %>% mutate_all(function(x) as.integer(x)) %>% as.matrix
}  

# write output for this chromosome  
saveRDS(sample_results, args$sample_results_out_path)
saveRDS(pair_results, args$pair_results_out_path)
cat("Saved results for chr", chr, ".\n")

# disconnect from CC database
dbDisconnect(snp.db)
