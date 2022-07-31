## code to prepare `internal` dataset goes here

lapply(list.files("data-raw/run_immucellai", pattern = "\\.rda$"), function(file) {
    load(paste0("data-raw/run_immucellai/", file), envir = globalenv())
})
immune_infiltrate_marker <- immune_infiltate_marker
compensation_matrix <- as.matrix(compensation_matrix)
typeof(compensation_matrix)
usethis::use_data(
    marker_exp, marker_exp_T, paper_marker,
    compensation_matrix, immune_infiltrate_marker,
    train_data, train_tag,
    internal = TRUE, overwrite = TRUE, 
    compress = "gzip"
)
