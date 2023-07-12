#' Prepare data for the input of PyClone or PyClone-vi
#'
#' @param mut_data A data.frame of mutation data.
#' @param cnv_data A data.frame of allele-specific CNV data with the top 5
#' columns containing "chromosome", "start_pos", "end_pos", "major_cn", and
#' "minor_cn". "sample_id" is optional, details see mut_data. other columns will
#' be omited. Column names don't matter.
#' @inheritParams identify_mut_cn
#' @param ref_field A string indicating the column names in `mut_data` that
#' contains the normal reference allele.
#' @param major_cn_field,minor_cn_field A string indicating the column names in
#' `cnv_data` that contains the major_cn and minor_cn.
#' @param ref_counts_field,var_counts_field A string indicating the column names
#' in `mut_data` that contains the ref_counts and var_counts.
#' @param normal_cn The copy number of the locus in non-malignant cells. This
#' should generally be 2 except for sex chromosomes in males.
#' @param pyclone_vi a logical value indicates whether prepare data for the
#' input of PyClone-vi. Details see <https://github.com/Roth-Lab/pyclone-vi>
#' @param purity_field Only used for
#' [PyClone-vi](https://github.com/Roth-Lab/pyclone-vi), a string specifying the
#' tumour content (cellularity) column (can in `mut_data` or `cnv_data`) of the
#' sample. Default value is 1.0 if NULL.
#' @param error_rate Only used for
#' [PyClone-vi](https://github.com/Roth-Lab/pyclone-vi), sequencing error rate.
#' @references
#'  - <https://github.com/Roth-Lab/pyclone>
#'  - <https://github.com/Roth-Lab/pyclone-vi>
#'  - <https://bitbucket.org/sequenzatools/sequenza/src/v2.1.1/R/next.R>
#' @export
prepare_pyclone <- function(
    mut_data, cnv_data, on_sample = NULL, on_chr = "chr",
    mut_pos = "pos", ref_field = "ref", start_field = "start",
    end_field = "end", major_cn_field = "major_cn", minor_cn_field = "minor_cn",
    ref_counts_field = "ref_counts",
    var_counts_field = "var_counts",
    purity_field = NULL, normal_cn = 2L,
    pyclone_vi = FALSE, error_rate = NULL,
    nomatch = NULL) {
    assert_df_with_columns(mut_data, c(
        names(on_sample) %||% on_sample,
        names(on_chr) %||% on_chr,
        mut_pos, ref_field,
        ref_counts_field, var_counts_field
    ))
    assert_df_with_columns(cnv_data, c(
        on_sample, on_chr, start_field, end_field,
        major_cn_field, minor_cn_field
    ))

    mut_data <- data.table::as.data.table(mut_data)
    cnv_data <- data.table::as.data.table(cnv_data)

    out <- mut_match_cn(
        mut_data, cnv_data,
        on_sample = on_sample,
        on_chr = on_chr, mut_pos = mut_pos,
        start_field = start_field, end_field = end_field,
        nomatch = nomatch
    )
    out[, mutation_id := Reduce(function(x, y) {
        paste(x, y, sep = ":")
    }, .SD), .SDcols = c(on_sample, on_chr, mut_pos, ref_field)]
    out[, normal_cn := normal_cn]
    columns <- c(
        "mutation_id", ref_counts_field, var_counts_field,
        "normal_cn", c(minor_cn_field, major_cn_field)
    )
    column_names <- c(
        "mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn",
        "major_cn"
    )
    if (pyclone_vi) {
        assert_class(purity_field, rlang::is_scalar_character,
            "scalar {.cls character}",
            null_ok = TRUE, cross_msg = NULL
        )
        if (is.null(purity_field)) {
            purity_field <- "..purity.."
            out$..purity.. <- 1L
        } else {
            assert_df_with_columns(out, purity_field,
                arg = c("mut_data", "cnv_data")
            )
        }
        columns <- c(columns, purity_field)
        column_names <- c(column_names, "purity")
    }

    out <- out[, .SD, SDcols = columns]
    data.table::setnames(out, column_names)
    if (pyclone_vi) {
        data.table::setnames(out, "var_counts", "alt_counts")
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
    if (!is.null(on_sample)) {
        data.table::setcolorder(out, "sample_id")
    }
    out[major_cn > 0L]
}
