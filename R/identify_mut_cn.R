#' Identify the copy number in mutation position
#' @param mut_data A data.frame of mutation data.
#' @param cnv_data A data.frame of allele-specific CNV data.
#' @param on_patient A string (can be named), specifying the patient column used
#' to match mut_data and cnv_data. If NULL, all mut_data and cnv_data will be
#' regarded from the same patient.
#' @param on_sample A string (which can be named), specifying the column in the
#' sample used for matching mut_data and cnv_data. If NULL, each patient is
#' considered to contain only one sample.
#' @param on_chr A string (can be named), specifying the chromosome column used
#' to match mut_data and cnv_data.
#' @param mut_pos A string indicating the column names in `mut_data` that
#' contains the variants positions in the chromosome.
#' @param start_field,end_field A string indicating the column names in
#' `cnv_data` that contains the start positions and end position of the genomic
#' ranges.
#' @param nomatch When a row in `mut_data` has no match to `cnv_data`,
#' nomatch=NA means NA is returned. NULL (or 0 for backward compatibility) means
#' no rows will be returned for that row of `mut_data`.
#' @return A integrated [data.table][data.table::data.table] with data column
#' from both mut_data and cnv_data.
#' @export
identify_mut_cn <- function(
    mut_data, cnv_data, on_patient = NULL, on_sample = NULL, on_chr = "chr",
    mut_pos = "pos", start_field = "start", end_field = "end", nomatch = NULL) {
    mut_match_cn(
        mut_data = mut_data, cnv_data = cnv_data,
        on_patient = on_patient, on_sample = on_sample,
        on_chr = on_chr, mut_pos = mut_pos,
        start_field = start_field, end_field = end_field,
        nomatch = nomatch
    )
}

#' @return A data.table
#' @keywords internal
#' @noRd
mut_match_cn <- function(
    mut_data, cnv_data, on_patient = NULL, on_sample = NULL, on_chr = "chr",
    mut_pos = "pos", start_field = "start", end_field = "end",
    nomatch = NULL, on_patient_arg = rlang::caller_arg(on_patient),
    on_sample_arg = rlang::caller_arg(on_sample),
    on_chr_arg = rlang::caller_arg(on_chr),
    mut_pos_arg = rlang::caller_arg(on_chr),
    start_field_arg = rlang::caller_arg(start_field),
    end_field_arg = rlang::caller_arg(end_field),
    call = parent.frame()) {
    assert_(on_patient, rlang::is_scalar_character,
        "a string",
        null_ok = TRUE, show_length = TRUE,
        arg = on_patient_arg,
        call = call
    )
    assert_(on_sample, rlang::is_scalar_character,
        "a string",
        cross_msg = NULL,
        null_ok = TRUE, show_length = TRUE,
        arg = on_sample_arg,
        call = call
    )
    assert_(on_chr, rlang::is_scalar_character,
        "a string",
        cross_msg = NULL, show_length = TRUE,
        arg = on_chr_arg,
        call = call
    )
    assert_(mut_pos, rlang::is_scalar_character,
        "a string",
        cross_msg = NULL, show_length = TRUE,
        arg = mut_pos_arg,
        call = call
    )
    assert_(start_field, rlang::is_scalar_character,
        "a string",
        cross_msg = NULL, show_length = TRUE,
        arg = start_field_arg,
        call = call
    )
    assert_(end_field, rlang::is_scalar_character,
        "a string",
        cross_msg = NULL, show_length = TRUE,
        arg = end_field_arg,
        call = call
    )
    mut_patient_col <- names(on_patient) %||% on_patient
    mut_sample_col <- names(on_sample) %||% on_sample
    mut_chr_col <- names(on_chr) %||% on_chr
    assert_data_frame(mut_data)
    assert_data_frame(cnv_data)
    assert_data_frame_columns(
        mut_data,
        c(mut_patient_col, mut_sample_col, mut_chr_col, mut_pos)
    )
    assert_data_frame_columns(
        cnv_data,
        c(on_sample, on_chr, start_field, end_field)
    )
    for (i in c(on_chr, start_field, end_field)) {
        if (anyNA(cnv_data[[i]])) {
            cli::cli_abort("{.val {NA}} is not allowed in {.code cnv_data[[{i}]]}")
        }
    }
    if (!is.null(on_patient)) {
        patients_in_mut_not_in_cnv <- setdiff(
            mut_data[[mut_patient_col]], cnv_data[[on_patient]]
        )
        patients_in_cnv_not_in_mut <- setdiff( # nolint
            cnv_data[[on_patient]], mut_data[[mut_patient_col]]
        )
        if (length(patients_in_mut_not_in_cnv)) {
            cli::cli_warn(c(
                "Cannot match all patients between {.arg mut_data} and {.arg cnv_data}",
                x = "Patients in {.arg mut_data} not in {.arg cnv_data}: {.val {patients_in_mut_not_in_cnv}}",
                i = "Patients in {.arg cnv_data} not in {.arg mut_data}: {.val {patients_in_cnv_not_in_mut}}"
            ))
        }
        on_string <- paste(on_patient, mut_patient_col, sep = "==")
    } else {
        on_string <- NULL
    }
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
        on_string <- c(on_string, paste(on_sample, mut_sample_col, sep = "=="))
    }

    on_string <- c(
        on_string, paste(on_chr, mut_chr_col, sep = "=="),
        paste("...start..___..pos...", "...mut..___..pos...", sep = "<="),
        paste("...end..___..pos...", "...mut..___..pos...", sep = ">=")
    )

    # Reduce the possibility of some columns in cnv_data be overrided
    # Also keep the original column mut_pos, start_field and end_field
    cnv_data <- data.table::as.data.table(cnv_data)
    cnv_data$...start..___..pos... <- cnv_data[[start_field]]
    # instead of use start_field directly since this column name may also exist
    # in mut_data, we just create a backup to know which mutation has not a
    # matched copy number, in this way, we can remove these rows.
    cnv_data$...start..___..pos...2 <- cnv_data[[start_field]]
    cnv_data$...end..___..pos... <- cnv_data[[end_field]]
    mut_data <- data.table::as.data.table(mut_data)
    mut_data$...mut..___..pos... <- mut_data[[mut_pos]]
    ntotal <- nrow(mut_data) # nolint
    out <- cnv_data[mut_data,
        on = on_string, nomatch = NA,
        allow.cartesian = FALSE
    ]

    # for every mutation, there must have one copy number value
    # abort if
    ntotal2 <- nrow(out)
    if (ntotal2 > ntotal) {
        msg <- "multiple copy number value found producing {ntotal2}/{ntotal} mutations"
        if (is.null(on_sample)) {
            msg <- c(msg, i = "try to set {.arg on_sample} or you should ensure no duplciated segment")
        } else {
            msg <- c(msg, i = "you should ensure no duplciated segment")
        }
        cli::cli_abort(msg)
    }
    # warning if nomatch found, since we have ensure start_field in cnv_data
    # have no NA value, it's save to regard NA as nomatch

    nomatch_rows <- is.na(out$...start..___..pos...2)
    if (any(nomatch_rows)) {
        cli::cli_warn(
            "Cannot match copy number for {sum(nomatch_rows)}/{ntotal} mutations"
        )
    }
    if (is.null(nomatch)) out <- out[!nomatch_rows]

    # remove the created columns
    out[, c("...start..___..pos...", "...start..___..pos...2", "...end..___..pos...") := list(NULL, NULL, NULL)]
    # check the match works well
    failed_pos <- out[[mut_pos]] < out[[start_field]] |
        out[[mut_pos]] > out[[end_field]]
    if (any(failed_pos, na.rm = TRUE)) {
        cli::cli_abort("Something wrong when parsing CN of mutation")
    }
    out[]
}
