# following code is modified based on
# <https://github.com/UCL-Research-Department-of-Pathology/panConusig/blob/main/R/setup_CNsigs.R>
#' Generate CN summary matrix from CN profiles
#' @inheritParams summarize_arm
#' @param minor_cn_field,major_cn_field A string specifying the minor_cn (minor
#'  allele) and major_cn (major allele) column in seg_data.
#' @references Christopher D. Steele, Maxime Tarabichi, et al, Undifferentiated
#' Sarcomas Develop through Distinct Evolutionary Pathways,Cancer Cell,
#' <https://doi.org/10.1016/j.ccell.2019.02.002>
#' @export
prepare_cnv_signature <- function(
    seg_data, sample_field = NULL, chr_field = "chr", start_field = "startpos",
    end_field = "endpos", major_cn_field = "major_cn", minor_cn_field = "minor_cn") {
    assert_df_with_columns(seg_data, c(
        sample_field, chr_field, start_field, end_field,
        minor_cn_field, major_cn_field
    ))
    seg_data <- data.table::as.data.table(seg_data)
    seg_data <- seg_data[, .SD, .SDcols = c(
        sample_field, chr_field, start_field, end_field,
        minor_cn_field, major_cn_field
    )]
    names <- c("chr", "startpos", "endpos", "minor_cn", "major_cn")
    if (!is.null(sample_field)) {
        names <- c("sample_id", names)
    }
    data.table::setnames(seg_data, names)

    # segment lengths
    width <- (seg_data$endpos - seg_data$startpos + 1L) / 1e6L # nolint
    CNt <- seg_data$minor_cn + seg_data$major_cn
    loh <- data.table::fcase(
        CNt == 0L, "homdel",
        pmin(seg_data$minor_cn, seg_data$major_cn) == 0L, "LOH",
        default = "het"
    )
    cn_class <- cut(CNt,
        breaks = c(-0.5, 0.5, 1.5, 2.5, 4.5, 8.5, Inf),
        labels = c("0", "1", "2", "3-4", "5-8", "9+")
    )
    seg_length_class <- data.table::fifelse(
        loh == "homdel",
        as.character(cut(width,
            breaks = c(-0.01, 0.1, 1, Inf)
        )),
        as.character(cut(width,
            breaks = c(-0.01, 0.1, 1, 10, 40, Inf)
        ))
    )
    length_to_class <- c(
        "(-0.01,0.1]" = "0-100kb",
        "(0.1,1]" = "100kb-1Mb",
        "(1,10]" = "1Mb-10Mb",
        "(10,40]" = "10Mb-40Mb",
        "(40,Inf]" = ">40Mb",
        "(1,Inf]" = ">1Mb"
    )
    seg_length_class <- length_to_class[seg_length_class]
    cnv_signatures <- paste(cn_class, loh, seg_length_class, sep = ":")
    table_counts(
        cnv_signatures, enumerate_standard_cnv_context(),
        sample_field %|n|% seg_data$sample_id
    )
}

enumerate_standard_cnv_context <- function() {
    out_list <- lapply(c("homdel", "LOH", "het"), function(loh) {
        cn_class <- switch(loh,
            homdel = "0",
            LOH = c("1", "2", "3-4", "5-8", "9+"),
            het = c("2", "3-4", "5-8", "9+")
        )
        length_class <- switch(loh,
            homdel = c("0-100kb", "100kb-1Mb", ">1Mb"),
            LOH = ,
            het = c("0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb")
        )
        cn_length_pair <- expand.grid(cn_class, length_class)
        paste(cn_length_pair[[1L]], loh, cn_length_pair[[2L]], sep = ":")
    })
    unlist(out_list, recursive = FALSE, use.names = FALSE)
}
