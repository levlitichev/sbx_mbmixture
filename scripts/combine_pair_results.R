pair_results_files <- snakemake@input
pair_results_out_path <- snakemake@output[[1]]

cat("Will combine", length(pair_results_files), "pair_results files.\n")

# initialize output
pair_results <- NULL

for (f in pair_results_files) {
    
    x <- readRDS(f)
    if(is.null(pair_results)) {
        pair_results <- x
    } else {
        if (length(x) != length(pair_results))
            stop(sprintf("file: %s, length(x): %i, length(pair_results): %i", f, length(x), length(pair_results)))
        stopifnot(dim(x) == dim(pair_results),
                  all(rownames(x) == rownames(pair_results)))
        pair_results <- pair_results + x
    }
}

# write output
saveRDS(pair_results, file=pair_results_out_path)
cat("Saved to", pair_results_out_path, "\n")
