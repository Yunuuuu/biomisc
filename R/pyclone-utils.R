#' Prepare data for the input of PyClone or PyClone-vi
#' @inheritParams identify_mut_cn
#' @param normal_cn The copy number of the locus in non-malignant cells. This
#' should generally be 2 except for sex chromosomes in males.
#' @param pyclone_vi a logical value indicates whether prepare data for the
#' input of PyClone-vi. Details see <https://github.com/Roth-Lab/pyclone-vi>
#' @param purity Only used for
#' [PyClone-vi](https://github.com/Roth-Lab/pyclone-vi), the tumour content
#' (cellularity) of the sample. Default value is 1.0 if column is not present.
#' @param error_rate Only used for
#' [PyClone-vi](https://github.com/Roth-Lab/pyclone-vi), sequencing error rate.
#' @references
#'  - <https://github.com/Roth-Lab/pyclone>
#'  - <https://github.com/Roth-Lab/pyclone-vi>
#'  - <https://bitbucket.org/sequenzatools/sequenza/src/v2.1.1/R/next.R>
#' @export
prepare_pyclone <- function(mut_data, cnv_data, sample_id = NULL, normal_cn = 2, pyclone_vi = FALSE, purity = 1L, error_rate = NULL) {
    if (pyclone_vi) {
        assert_class(purity, is_scalar_numeric,
            "scalar numeric",
            cross_msg = NULL
        )
    }
    out <- identify_mut_cn(mut_data, cnv_data, sample_id = sample_id)
    out[, mutation_id := paste(chromosome, position, sep = ":")]
    out[, normal_cn := normal_cn]

    data.table::setcolorder(
        out,
        c(
            "mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn",
            "major_cn"
        )
    )
    if (pyclone_vi) {
        # out[, sample_id := sample_id]
        out[, tumour_content := purity]
        data.table::setnames(out, "var_counts", "alt_counts")
        # data.table::setcolorder(out,
        #     "sample_id",
        #     after = "mutation_id"
        # )
        data.table::setcolorder(out,
            c("minor_cn", "normal_cn"),
            after = "major_cn"
        )
        data.table::setcolorder(out,
            "tumour_content",
            after = "normal_cn"
        )
        if (!is.null(error_rate)) {
            out[, error_rate := error_rate]
            data.table::setcolorder(out,
                "error_rate",
                after = "tumour_content"
            )
        }
    }
    if (!is.null(sample_id)) data.table::setcolorder(out, "sample_id")
    out[major_cn > 0L]
}

#' Identify the copy number in mutation position
#' @param mut_data A data.frame of mutation data with the top 4 columns (omit
#' column specified in `sample_id`) containing "chromosome", "pos", "ref_counts"
#' and "var_counts", other columns will be omited. Column names don't matter. If
#' `sample_id` is not NULL, the column specified in `sample_id` will be matched
#' in both mut_data and cnv_data.
#' @param cnv_data A data.frame of allele-specific CNV data with the top 5
#' columns containing "chromosome", "start_pos", "end_pos", "major_cn", and
#' "minor_cn". "sample_id" is optional, details see mut_data. other columns will
#' be omited. Column names don't matter.
#' @param sample_id A string (can be named), specifying the column used to match
#' mut_data and cnv_data.
#' @export
identify_mut_cn <- function(mut_data, cnv_data, sample_id = NULL) {
    assert_class(mut_data, function(x) {
        inherits(x, "data.frame") && ncol(x) >= 4L
    }, "{.cls data.frame} with at least 4 columns")
    assert_class(cnv_data, function(x) {
        inherits(x, "data.frame") && ncol(x) >= 5L
    }, "{.cls data.frame} with at least 5 columns")
    assert_class(sample_id, rlang::is_scalar_character,
        "scalar character",
        null_ok = TRUE
    )
    if (!is.null(sample_id)) {
        cnv_sample_col <- unname(sample_id)
        mut_sample_col <- names(sample_id) %||% cnv_sample_col
        cnv_sample_id <- cnv_data[[cnv_sample_col]]
        cnv_data[[cnv_sample_col]] <- NULL
        mut_sample_id <- mut_data[[mut_sample_col]]
        mut_data[[mut_sample_col]] <- NULL
    } else {
        cnv_sample_id <- rep_len("sample", nrow(cnv_data))
        mut_sample_id <- rep_len("sample", nrow(mut_data))
    }

    mut_data <- data.table::as.data.table(mut_data)[, 1:4]
    data.table::setnames(
        mut_data,
        c("chromosome", "pos", "ref_counts", "var_counts")
    )
    mut_data[, sample_id := mut_sample_id]

    cnv_data <- data.table::as.data.table(cnv_data)[, 1:5]
    data.table::setnames(
        cnv_data,
        c("chromosome", "start_pos", "end_pos", "major_cn", "minor_cn")
    )
    cnv_data[, sample_id := cnv_sample_id]

    mut_cn <- cnv_data[
        mut_data,
        list(
            sample_id = i.sample_id,
            ref_counts = i.ref_counts, var_counts = i.var_counts,
            minor_cn = x.minor_cn, major_cn = x.major_cn,
            chromosome = i.chromosome, position = i.pos,
            start_pos = x.start_pos, end_pos = x.end_pos
        ),
        on = c("sample_id", "chromosome", "start_pos<=pos", "end_pos>=pos"),
        nomatch = NULL
    ][, .SD[rowSums(is.na(.SD)) == 0L]]
    failed_pos <- mut_cn[["pos"]] < mut_cn[["start_pos"]] |
        mut_cn[["pos"]] > mut_cn[["end_pos"]]
    if (any(failed_pos)) {
        cli::cli_abort("Something wrong when parsing CN of mutation")
    }
    if (is.null(sample_id)) mut_cn$sample_id <- NULL
    mut_cn
}
