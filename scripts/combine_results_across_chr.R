sample_results_files <- snakemake@input[["sample_results"]]
pair_results_files <- snakemake@input[["pair_results"]]
sample_results_out_path <- snakemake@output[["sample_results"]]
pair_results_out_path <- snakemake@output[["pair_results"]]

cat("Will combine", length(sample_results_files), "sample_results files.\n")
cat("Will combine", length(pair_results_files), "pair_results files.\n")

# initialize outputs
pair_results <- sample_results <- NULL

# sample results
for (f in sample_results_files) { 
    
    x <- readRDS(f)
    if(is.null(sample_results)) {
        sample_results <- x
    } else {
        if (length(x) != length(sample_results))
            stop(sprintf("file: %s, length(x): %i, length(sample_results): %i", f, length(x), length(sample_results)))
        stopifnot(dim(x) == dim(sample_results),
                  all(rownames(x) == rownames(sample_results)))
        sample_results <- sample_results + x
    }
}

# write output
saveRDS(sample_results, file=sample_results_out_path)
cat("Saved to", sample_results_out_path, "\n")

# pair results
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
