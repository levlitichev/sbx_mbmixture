library(mbmixture) # github.com/kbroman/mbmixture

# this script is currently run manually, not as part of snakemake
sample.results.files <- Sys.glob("/path/to/sunbeam/sunbeam_output/qc/mbmixture/genome/summaries/*/sample_results_allchr.rds")
pair.results.files <- Sys.glob("/path/to/sunbeam/sunbeam_output/qc/mbmixture/genome/summaries/*/pair_results_allchr.rds")
sample.results.out.path <- "/path/to/sunbeam/sunbeam_output/qc/mbmixture/genome/summaries/sample_results_all.csv"
pair.results.out.path <- "/path/to/sunbeam/sunbeam_output/qc/mbmixture/genome/summaries/pair_results_all.csv"

# sample_results
list.of.sample.result.summaries <- vector("list", length=length(sample.results.files))
for (ii in 1:length(sample.results.files)) {
    this.res <- readRDS(sample.results.files[[ii]])
    this.mb.sample <- basename(dirname(sample.results.files[[ii]]))
    
    # calculate proportion of mismatches at homozygous loci
    mismatches <- apply(this.res, 1, function(mat) {mat[1,2] + mat[3,1]})
    total.counts.at.homo.loci <- apply(this.res, 1, function(mat) {sum(mat[1,] + mat[3,])})
    this.sample.result.summary <- data.frame(
        mb.sample=this.mb.sample,
        num.mismatches=mismatches,
        total.counts.at.homo.loci=total.counts.at.homo.loci)
    this.sample.result.summary$mouse.ID <- rownames(this.sample.result.summary)
    rownames(this.sample.result.summary) <- NULL

    list.of.sample.result.summaries[[ii]] <- this.sample.result.summary 
}

sample.result.summary.df <- data.frame(do.call(rbind, list.of.sample.result.summaries))

# write output
write.csv(apply(sample.result.summary.df, 2, as.character),
          file=sample.results.out.path, quote=F, row.names=F)
cat("Wrote", sample.results.out.path, "\n")

# pair_results
list.of.pair.result.summaries <- vector("list", length=length(pair.results.files))
for (ii in 1:length(pair.results.files)) {
    
    this.pair.res <- readRDS(pair.results.files[[ii]])
    cat("Loaded", pair.results.files[[ii]], "\n")
    this.mb.sample <- basename(dirname(pair.results.files[[ii]]))

    this.pair.summary.df <- data.frame(t(apply(this.pair.res, 1, mle_pe)))
    this.pair.summary.df$mouse.ID <- rownames(this.pair.summary.df)
    this.pair.summary.df$mb.sample <- this.mb.sample
    rownames(this.pair.summary.df) <- NULL

    list.of.pair.result.summaries[[ii]] <- this.pair.summary.df
}

pair.result.summary.df <- do.call(rbind, list.of.pair.result.summaries)

# write output
write.csv(apply(pair.result.summary.df, 2, as.character),
          pair.results.out.path,  quote=F, row.names=F)
cat("Wrote", pair.results.out.path, "\n")
