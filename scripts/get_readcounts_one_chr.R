suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(optparse))

# check if script is being run via snakemake
if (!exists("snakemake"))
    stop(sprintf("This script only works as part of a Snakemake pipeline.")) 
    
bam_path <- snakemake@input[[1]]
out_path <- snakemake@output[[1]]
snp_db_path <- snakemake@config[["sbx_mbmixture"]][["CC_sqlite_db"]]
chr <- as.integer(snakemake@params[["chr"]])
bam_has_chr_prefix <- TRUE #TODO

# make sure bam file exists
if (!file.exists(bam_path))
    stop(sprintf("Path to bam file doesn't exist. bam_path: %s", bam_path))

# make sure index file exists next to bam
if (!file.exists(paste0(bam_path, ".bai")))
    stop(sprintf("Can't find index file at %s", paste0(bam_path, ".bai"))) 

# make sure chromosome was provided
if (is.null(chr))
    stop("Must provide a chromosome.")
cat("Chromosome", chr, "\n")

# make sure CC SNP db exists
if (!file.exists(snp_db_path))
    stop(sprintf("Path to CC SNP db doesn't exist. snp_db_path: %s", snp_db_path))

# make sure output path can be written to
if (!dir.exists(dirname(out_path)))
    stop(sprintf("Can't write to the output path. out_path: %s", out_path))

# load bam file (index must be right next to the file)
bf <- BamFile(bam_path, index=paste0(bam_path, ".bai"))
bf.seqinfo <- as.data.frame(seqinfo(bf))
cat("Loaded", bam_path, "\n")

# connect to CC database
snp.db <- dbConnect(SQLite(), snp_db_path)

# optionally add "chr" prefix
if (bam_has_chr_prefix) {chr.for.bam.file <- paste0("chr", chr)
} else {chr.for.bam.file <- chr}

# compute pileup at each location (for this chromosome in this bam file)
chr_length <- bf.seqinfo[chr.for.bam.file, "seqlengths"]
scanBamParam <- ScanBamParam(which=setNames( IRangesList(IRanges(0L, chr_length)), chr.for.bam.file) )
pileup.res <- pileup(bf, scanBamParam=scanBamParam) # slow

# extract positions and possible alleles for CC SNPs 
CC_snp_info <- dbGetQuery(snp.db, paste0("select pos,alleles from variants where chr=='", chr, "'"))
tmp <- strsplit(CC_snp_info$alleles, "\\|")
CC_snp_info$allele1 <- sapply(tmp, "[", 1)
CC_snp_info$allele2 <- sapply(tmp, "[", 2)

# only keep biallelic SNPs
CC_snp_info <- CC_snp_info[
  CC_snp_info$allele1 %in% c("A", "C", "G", "T") & CC_snp_info$allele2 %in% c("A", "C", "G", "T"), ]

# only keep pileup results for bi-allelic CC SNPs
pileup.res <- pileup.res[pileup.res$pos %in% CC_snp_info$pos,,drop=FALSE]

stopifnot(nrow(pileup.res) > 0)
cat("Found", sum(pileup.res$count), "read counts at", length(unique(pileup.res$pos)), "bi-allelic CC SNPs\n")

# initialize output df
out.counts.df <- CC_snp_info[CC_snp_info$pos %in% unique(pileup.res$pos), ]
out.counts.df$count2 <- out.counts.df$count1 <- 0
out.counts.df$chr <- chr

# for every bi-allelic CC SNP, count how many reads were found for each of the two alleles 
tmp.counts.df <- tapply(pileup.res$count, list(pileup.res$pos, pileup.res$nucleotide), sum)
tmp.counts.df[is.na(tmp.counts.df)] <- 0
m <- match(out.counts.df$pos, rownames(tmp.counts.df))
for(a in c("A", "C", "G", "T")) {
    out.counts.df[out.counts.df$allele1==a,"count1"] <- tmp.counts.df[m,][out.counts.df$allele1==a,a]
    out.counts.df[out.counts.df$allele2==a,"count2"] <- tmp.counts.df[m,][out.counts.df$allele2==a,a]
}

# disconnect from CC database
dbDisconnect(snp.db)

# write output
out.df <- out.counts.df[, c("chr", "pos", "count1", "count2")]
write.csv(out.df, gzfile(out_path), quote=F, row.names=F)
cat("Output written to", out_path, "\n")
