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

utils::data(diseaseMap, # nolint
    package = "ABSOLUTE", envir = environment()
)
absolute_disease_map_env <- get(
    "disease_map",
    envir = environment()
)
absolute_disease_map_data <- ls(
    envir = absolute_disease_map_env, all.names = TRUE
)
cellmarker_database <- list(
    human = data.table::fread(
        "http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt"
    ),
    mouse = data.table::fread(
        "http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt"
    )
)
cellmarker_database <- lapply(cellmarker_database, cellmarker_prepare)

usethis::use_data(
    cibersort_lm22,
    absolute_disease_map_data,
    cellmarker_database,
    internal = TRUE, overwrite = TRUE,
    compress = "gzip"
)
