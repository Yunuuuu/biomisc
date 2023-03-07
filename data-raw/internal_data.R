## code to prepare `internal` dataset goes here

# run_arm_cnv --------------------------------------------------------

# ** ref_cytoband ---------------------------------------------------------

anno_hub <- AnnotationHub::AnnotationHub(localHub = TRUE)
AnnotationHub::query(
    anno_hub, c("UCSC", "Homo sapiens", "cytoband"),
    ignore.case = TRUE
)
# AH53177 | UCSC cytoBand track for hg19
# AH53178 | UCSC cytoBand track for hg38

run_arm_cnv_ref_cytoband_hg19 <- anno_hub[["AH53177"]]
run_arm_cnv_ref_cytoband_hg38 <- anno_hub[["AH53178"]]

# run_cibersort -----------------------------------------------------------
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

usethis::use_data(
    cibersort_lm22,
    absolute_disease_map_data,
    run_arm_cnv_ref_cytoband_hg19,
    run_arm_cnv_ref_cytoband_hg38,
    internal = TRUE, overwrite = TRUE,
    compress = "gzip"
)
