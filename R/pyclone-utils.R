#' Prepare data for the input of PyClone or PyClone-vi
#'
#' @param mut_data A data.frame of mutation data.
#' @param cnv_data A data.frame of allele-specific CNV data with the top 5
#' columns containing "chromosome", "start_pos", "end_pos", "major_cn", and
#' "minor_cn". "sample_id" is optional, details see mut_data. other columns will
#' be omited. Column names don't matter.
#' @inheritParams identify_mut_cn
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
    mut_pos = "pos", start_field = "start", end_field = "end",
    major_cn_field = "major_cn", minor_cn_field = "minor_cn",
    ref_counts_field = "ref_counts", var_counts_field = "var_counts",
    purity_field = NULL, normal_cn = 2L,
    pyclone_vi = FALSE, error_rate = NULL,
    nomatch = NA) {
    assert_df_with_columns(mut_data, c(
        names(on_sample) %||% on_sample,
        names(on_chr) %||% on_chr,
        mut_pos, ref_counts_field, var_counts_field
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
    out[, mutation_id := paste(chromosome, position, sep = ":")]
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
            null_ok = TRUE
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

#' Identify the copy number in mutation position
#' @param mut_data A data.frame of mutation data.
#' @param cnv_data A data.frame of allele-specific CNV data.
#' @param on_sample A string (can be named), specifying the sample column used
#' to match mut_data and cnv_data. If NULL, all mut_data and cnv_data will be
#' regarded from the same sample.
#' @param on_chr A string (can be named), specifying the chromosome column used
#' to match mut_data and cnv_data.
#' @param mut_pos A string indicating the column names in `mut_data` that
#' contains the variants positions in the chromosome.
#' @param start_field,end_field A string indicating the column names in
#' `cnv_data` that contains the start positions and end position of the genomic
#' ranges.
#' @param nomatch When a row in `mut_data` has no match to `cnv_data`,
#' nomatch=NA (default) means NA is returned. NULL (or 0 for backward
#' compatibility) means no rows will be returned for that row of `mut_data`.
#' @return A integrated data.frame with data column from both mut_data and
#' cnv_data.
#' @export
identify_mut_cn <- function(
    mut_data, cnv_data, on_sample = NULL, on_chr = "chr", mut_pos = "pos",
    start_field = "start", end_field = "end", nomatch = NA) {
    mut_cn <- mut_match_cn(
        mut_data = mut_data, cnv_data = cnv_data,
        on_sample = on_sample, on_chr = on_chr,
        mut_pos = mut_pos, start_field = start_field,
        end_field = end_field, nomatch = nomatch
    )
    data.table::setDF(mut_cn)
    mut_cn
}

#' @return A data.table
#' @keywords internal
#' @noRd
mut_match_cn <- function(
    mut_data, cnv_data, on_sample = NULL, on_chr = "chr",
    mut_pos = "pos", start_field = "start", end_field = "end",
    nomatch = NA,
    on_sample_arg = rlang::caller_arg(on_sample),
    on_chr_arg = rlang::caller_arg(on_chr),
    mut_pos_arg = rlang::caller_arg(on_chr),
    start_field_arg = rlang::caller_arg(start_field),
    end_field_arg = rlang::caller_arg(end_field),
    call = parent.frame()) {
    assert_class(on_sample, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE,
        arg = on_sample_arg,
        call = call
    )
    assert_class(on_chr, rlang::is_scalar_character,
        "scalar {.cls character}",
        arg = on_chr_arg,
        call = call
    )
    assert_class(mut_pos, rlang::is_scalar_character,
        "scalar {.cls character}",
        arg = mut_pos_arg,
        call = call
    )
    assert_class(start_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        arg = start_field_arg,
        call = call
    )
    assert_class(end_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        arg = end_field_arg,
        call = call
    )
    mut_sample_col <- names(on_sample) %||% on_sample
    mut_chr_col <- names(on_chr) %||% on_chr
    assert_df_with_columns(mut_data, c(mut_sample_col, mut_chr_col, mut_pos))
    assert_df_with_columns(
        cnv_data,
        c(on_sample, on_chr, start_field, end_field)
    )
    if (!is.null(on_sample)) {
        samples_in_mut_not_in_cnv <- setdiff(
            mut_data[[mut_sample_col]], cnv_data[[on_sample]]
        )
        samples_in_cnv_not_in_mut <- setdiff( # nolint
            cnv_data[[on_sample]], mut_data[[mut_sample_col]]
        )
        if (length(samples_in_mut_not_in_cnv)) {
            cli::cli_warn(c(
                "Cannot match all samples between {.arg mut_data} and {.arg cnv_data}",
                x = "Samples in {.arg mut_data} not in {.arg cnv_data}: {.val {samples_in_mut_not_in_cnv}}",
                i = "Samples in {.arg cnv_data} not in {.arg mut_data}: {.val {samples_in_cnv_not_in_mut}}"
            ))
        }
        on_string <- paste(on_sample, mut_sample_col, sep = "==")
    } else {
        on_string <- character()
    }
    # Reduce the possibility of some columns in cnv_data named as mut_data
    on_string <- c(
        on_string, paste(on_chr, mut_chr_col, sep = "=="),
        paste(start_field, "...mut_pos...", sep = "<="),
        paste(end_field, "...mut_pos...", sep = ">=")
    )
    ..mut_data.. <- data.table::as.data.table(mut_data)
    ..mut_data..$...mut_pos... <- ..mut_data..[[mut_pos]]
    mut_cn <- data.table::as.data.table(cnv_data)[
        ..mut_data..,
        on = on_string, nomatch = nomatch
    ]
    # check the match works well
    failed_pos <- mut_cn[[mut_pos]] < mut_cn[[start_field]] |
        mut_cn[[mut_pos]] > mut_cn[[end_field]]
    if (any(failed_pos)) {
        cli::cli_abort("Something wrong when parsing CN of mutation")
    }
    mut_cn
}
