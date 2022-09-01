## code to prepare `internal` dataset goes here

cibersort_lm22 <- readr::read_tsv(
    "data-raw/run_cibersort/LM22.txt",
    col_names = TRUE
)
cibersort_lm22 <- data.table::setDF(
    cibersort_lm22[-1L], 
    rownames = cibersort_lm22[["Gene symbol"]]
)
cibersort_lm22 <- as.matrix(cibersort_lm22)
usethis::use_data(
    cibersort_lm22,
    internal = TRUE, overwrite = TRUE, 
    compress = "gzip"
)
