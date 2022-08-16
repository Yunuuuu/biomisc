## code to prepare `internal` dataset goes here

usethis::use_data(
    marker_exp, marker_exp_T, paper_marker,
    compensation_matrix, immune_infiltrate_marker,
    train_data, train_tag,
    internal = TRUE, overwrite = TRUE, 
    compress = "gzip"
)
