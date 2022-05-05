sample_results_files <- snakemake@input[[1]]
sample_results_out_path <- snakemake@output[[1]]

cat("Will combine", length(sample_results_files), "sample_results files.\n")

# initialize output
sample_results <- NULL

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
