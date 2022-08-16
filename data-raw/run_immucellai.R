if (!dir.exists("inst/extdata/immucellai ")) {
    dir.create("inst/extdata/immucellai ", recursive = TRUE)
}

sample <- read.table(
    "data-raw/run_immucellai/example.txt",
    header = TRUE
)
sample <- sample[-1L, , drop = FALSE]
rownames(sample) <- sample[, 1L, drop = TRUE]
sample <- sample[, -1L, drop = FALSE]
sample_exp <- apply(as.matrix(sample), 2, as.numeric)
rownames(sample_exp) <- rownames(sample)

saveRDS(
    sample_exp, "inst/extdata/immucellai/sample_exp.rds",
    compress = "gzip"
)

lapply(list.files("data-raw/run_immucellai", pattern = "\\.rda$"), function(file) {
    load(paste0("data-raw/run_immucellai/", file), envir = globalenv())
    if (!file %in% c("immune_infiltate_marker.rda", "compensation_matrix.rda")) {
        saveRDS(
            rlang::env_get(
                globalenv(),
                nm = fs::path_ext_remove(basename(file))
            ),
            paste0(
                "inst/extdata/immucellai/",
                fs::path_ext_remove(basename(file)), ".rds"
            ),
            compress = "gzip"
        )
    }
})
immune_infiltrate_marker <- immune_infiltate_marker
saveRDS(
    immune_infiltrate_marker,
    paste0("inst/extdata/immucellai/", "immune_infiltrate_marker", ".rds"),
    compress = "gzip"
)
compensation_matrix <- as.matrix(compensation_matrix)
saveRDS(
    compensation_matrix,
    paste0("inst/extdata/immucellai/", "compensation_matrix", ".rds"),
    compress = "gzip"
)
