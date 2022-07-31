sample <- read.table(
    "data-raw/run_immucellai/example.txt",
    header = TRUE
)
sample <- sample[-1L, , drop = FALSE]
rownames(sample) <- sample[, 1L, drop = TRUE]
sample <- sample[, -1L, drop = FALSE]
sample_exp <- apply(as.matrix(sample), 2, as.numeric)
rownames(sample_exp) <- rownames(sample)

if (!dir.exists("inst/extdata")) {
    dir.create("inst/extdata", recursive = TRUE)
}
saveRDS(
    sample_exp, "inst/extdata/run_immucellai_sample_exp.rds",
    compress = "gzip"
)

