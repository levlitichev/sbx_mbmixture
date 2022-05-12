library(mbmixture) # github.com/kbroman/mbmixture

# this script is currently run manually, not as part of snakemake
pair.results.files <- Sys.glob("/path/to/sunbeam/sunbeam_output/qc/mbmixture/mm10/summaries/pair_results/*/all.rds")
pair.results.out.path <- "/path/to/sunbeam/sunbeam_output/qc/mbmixture/mm10/summaries/all_pair_results.csv"

list.of.pair.result.summaries <- vector("list", length=length(pair.results.files))
for (ii in 1:length(pair.results.files)) {
    
    this.pair.res <- readRDS(pair.results.files[[ii]])
    cat("Loaded", pair.results.files[[ii]], "\n")
    this.mb.sample <- basename(dirname(pair.results.files[[ii]]))

    this.pair.summary.df <- data.frame(t(apply(this.pair.res, 1, mle_pe))) # mle_pe from mbmixture
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
