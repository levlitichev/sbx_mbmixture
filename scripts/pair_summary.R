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

##### WARNING #####
# The following line of code is fragile and depends on how mouse ID is encoded in the sample name

# if the sample is named like `DO_1D_3011_044w` and the mouse should be `DO-1D-3011`:
#mouse.ID.for.this.mb.sample <- paste0(strsplit(sample, "_")[[1]][1:3], collapse="-")

# if the sample is named like `LB17_DO_1D_3011_044w` and the mouse should be `DO-1D-3011`:
mouse.ID.for.this.mb.sample <- paste0(strsplit(sample, "_")[[1]][2:4], collapse="-")

cat("Mouse ID for this microbiome sample:", mouse.ID.for.this.mb.sample, "\n")
##### END OF WARNING #####

# read just one line of imp_snps to get mouse_IDs
imp_snps_tmp <- read.csv(imp_snp_path, nrow=1, header=T, check.names=F)
mouse_IDs <- colnames(imp_snps_tmp)[3:ncol(imp_snps_tmp)]

# check that the mouse ID for our microbiome sample matches a mouse ID in imp_snps
if (!(mouse.ID.for.this.mb.sample %in% mouse_IDs))
  stop(sprintf("No mouse ID in imp_snps matches our mouse ID. Mouse ID in microbiome sample: %s, first 3 mice in imp_snps: %s", mouse.ID.for.this.mb.sample, paste0(mouse_IDs[1:3], collapse=", ")))

# import all of imp_snps (slow)
imp_snps <- read_csv(imp_snp_path, col_types=cols(chr="i", pos="i", .default="d"), num_threads=threads)
cat("Loaded", imp_snp_path, "\n")

# SNPs in counts.df should be a subset of imputed bi-allelic SNPs
# (because in get_readcounts.R, we already subset to bi-allelic SNPs in the CC sqlite db)
merged  <- merge(imp_snps, counts.df, by=c("chr", "pos"))
stopifnot(nrow(merged) > 0)
cat("Using", nrow(merged), "SNPs to compare this microbiome sample to genotyped mice...\n")

# pair_results is used for determining whether a microbiome sample is actually a mixture of two samples
# for all genotypes, simultaneously gets counts for the expected mouse genotype versus another genotype

# initialize output
pair_results <- array(0, dim=c(length(mouse_IDs), 3, 3, 2))
dimnames(pair_results) <- list(mouse_IDs, c("AA", "AB", "BB"),
                               c("AA", "AB", "BB"), c("A", "B"))
# loop over mice
for (ii in 1:length(mouse_IDs)) {
  cat("Comparing microbiome sample to mouse", ii, "of", length(mouse_IDs), "\n") 
  this.mouse <- mouse_IDs[ii]
  
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

# write output  
saveRDS(pair_results, out_path)
cat("Saved pair_results for chr", chr, "\n")
