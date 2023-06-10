#' Matching mutation Copy number and Estimating CCF
#'
#' @description
#' Following methods to define subclone and clone.
#' - define clone and subclone based on absCCF or mut_multi_btstr
#' - define early or late with Mt or phyloCCF.
#'
#' For CONIPHER anlayis, use min_subclonal = 0.05, conipher = TRUE,
#' min_vaf_to_explain = 0.05. And you should provide indel_field or both
#' ref_field and alt_field.
#'
#' @inheritParams identify_mut_cn
#' @param ccf_type which CCF should we estimate, one of "phyloCCF" (include both
#' absolute and phylogenetic CCF) or "bootstrapCCF", or "both".
#' @param purity_field A string specifying the purity column (in `mut_data` or
#'  `cnv_data`). Default is "purity".
#' @param normal_cn A scalar number specifying the normal.copy number or a
#'  string indicating the normal_cn column in `mut_data` since it define normal
#'  copy number for every mutation. It's save to use 2 if you only analyze
#'  autosomes. Or you should use `gender_field` to define the `normal_cn`. For
#'  sex chromosomes and gender is "male", 1L will be used, otherwise, 2L will be
#'  used.
#' @param gender_field A string specifying the gender column. Only used when
#'   normal_cn is NULL. Default is "gender". Only "female" and "male" are
#'   supported in this column.
#' @param contigs An atomic vector specifying the chromosome to analyze. If
#' NULL, all chromosome will be used.
#' @param ... Not used currently
#' @param ref_field,alt_field A string specifying the column of reference allele
#'  or variant allele in mut_data to estimate whether a variant is a indel or
#'  not. `alt_field` will only be used when `indel_field` is NULL. `ref_field`
#'  must be specified when `indel_field` or `alt_field` is not NULL.
#' @param indel_field A string specifying a logical column of mut_data which
#' indicates whether or not a variant is a indel.
#' @param min_vaf_to_explain A numeric, the minimal vaf value to define
#'  subclone.
#' @param min_subclonal Minimal copy number to define subclone.
#' @param conipher A scalar logical indicates whether calculate phyloCCF like
#'  CONIPHER. Details see
#'  <https://github.com/McGranahanLab/CONIPHER-wrapper/tree/main>
#' @param nomatch When a row in `mut_data` has no match to `cnv_data`,
#' nomatch=NA means NA is returned. NULL (or 0 for backward compatibility) means
#' no rows will be returned for that row of `mut_data`.  It's not unusual to
#' place purity data in cnv_data, in this way, if some mutation cannot match the
#' cnv data, the purity and copy number value for this mutation would be NA, For
#' CCF estimation, NA is not allowed. Just set nomatch = NULL to omit these
#' rows.
#' @param kept_cols A character vector specifying the columns in `mut_data` or
#' `cnv_data` you want to keep in the results. By default only used column and
#' created columns will be returned.
#' @seealso
#'  - <https://bitbucket.org/nmcgranahan/clonalneoantigenanalysispipeline>
#'  - <https://bitbucket.org/nmcgranahan/pancancerclonality>
#'  - <https://github.com/McGranahanLab/CONIPHER-wrapper/>
#' @note
#' `on_sample` argument is just for multi-region from a single patient. You
#' should always use `on_patient` to specify the matched sample for a
#' single-region protocol if you use phyloCCF and want to correct indel
#' phyloCCF.
#' @export
run_ccf <- function(
    mut_data, cnv_data,
    ccf_type = c("phyloCCF", "bootstrapCCF", "both"),
    on_patient = NULL, on_sample = NULL,
    purity_field = NULL, on_chr = "chr", mut_pos = "pos",
    start_field = "startpos", end_field = "endpos",
    normal_cn = 2L, gender_field = NULL, contigs = 1:22,
    ...,
    # arguments for ccf_type = phyloCCF
    ref_field = NULL,
    alt_field = NULL, indel_field = NULL,
    min_vaf_to_explain = NULL, min_subclonal = NULL, conipher = FALSE,
    # Other arguments
    nomatch = NULL, kept_cols = NULL) {
    ccf_type <- match.arg(ccf_type)
    rlang::check_dots_empty()

    # check arguments firstly
    assert_class(on_patient, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE, cross_msg = NULL
    )
    assert_class(on_sample, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE, cross_msg = NULL
    )
    assert_class(purity_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE, cross_msg = NULL
    )
    assert_class(
        on_chr, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    assert_class(
        mut_pos, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    assert_class(
        start_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    assert_class(
        end_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    assert_class(normal_cn,
        function(x) {
            is_scalar_numeric(x) || rlang::is_scalar_character(x)
        }, "scalar {.cls numeric} or {.cls character}",
        null_ok = TRUE, cross_msg = NULL
    )

    assert_df_with_columns(mut_data, c(
        names(on_patient) %||% on_patient,
        names(on_sample) %||% on_sample,
        names(on_chr) %||% on_chr, mut_pos,
        "ref_counts", "alt_counts"
    ))
    cnv_columns <- switch(ccf_type,
        phyloCCF = c("nAraw", "nBraw"),
        bootstrapCCF = c("nMinor", "nMajor"),
        both = c("nAraw", "nBraw", "nMinor", "nMajor")
    )
    assert_df_with_columns(cnv_data, c(
        on_patient, on_sample, on_chr, start_field, end_field, cnv_columns
    ))
    assert_class(kept_cols, is.character, "{.cls character} vector",
        cross_msg = NULL, null_ok = TRUE
    )
    mut_data <- data.table::as.data.table(mut_data)
    cnv_data <- data.table::as.data.table(cnv_data)

    if (ccf_type != "bootstrapCCF") {
        # define subclone copy number
        cnv_data <- define_subclone_cn(cnv_data, min_subclonal = 0.01)
    }

    # just extract the segmented CNV for this sample
    out <- mut_match_cn(mut_data, cnv_data,
        on_patient = on_patient, on_sample = on_sample, on_chr = on_chr,
        mut_pos = mut_pos, start_field = start_field,
        end_field = end_field, nomatch = nomatch
    )

    # assert every sample has only one patient
    if (!is.null(on_sample)) assert_nest(out, on_patient, on_sample)
    if (!is.null(on_patient) && is.null(on_sample)) {
        group <- on_patient
    } else {
        group <- on_sample
    }
    if (is.null(group)) {
        info_msg <- "try to set {.arg on_patient} or {.arg on_sample}"
    } else {
        info_msg <- NULL
    }

    # assert every necessary column is in the combined data
    purity_field <- purity_field %||% "purity"
    assert_df_with_columns(out, c(purity_field, kept_cols),
        check_class = FALSE, arg = c("mut_data", "cnv_data")
    )

    # assert every samples provided only have one purity value
    assert_nest(out, purity_field, group, info_msg = info_msg)
    data.table::setnames(out, purity_field, "purity")
    if (!all(data.table::between(out$purity, 0L, 1L))) {
        cli::cli_abort("purity must in [0, 1]")
    }
    if (!is.null(contigs)) {
        # filter contigs
        matched_contigs <- as.character(out[[on_chr]]) %chin%
            as.character(contigs)
        out <- out[matched_contigs]
    }

    # define normal_cn
    if (is.null(normal_cn)) {
        # assert every samples provided only one gender value
        gender_field <- gender_field %||% "gender"
        assert_df_with_columns(out, gender_field,
            check_class = FALSE,
            arg = c("mut_data", "cnv_data")
        )
        if (!all(out[[gender_field]] %in% c("male", "female"))) {
            cli::cli_abort("Only {.val male} and {.val female} are supported in {.field {gender_field}} column")
        }
        assert_nest(
            out, gender_field, group,
            cross_format = "group", info_msg = info_msg
        )
        out$normal_cn <- define_normal_cn(out[[gender_field]], out[[on_chr]])
    } else {
        if (is_scalar_numeric(normal_cn)) {
            out[, normal_cn := normal_cn]
        } else {
            assert_df_with_columns(out, normal_cn,
                check_class = FALSE,
                arg = "mut_data"
            )
            if (!is.numeric(out[[normal_cn]])) {
                cli::cli_abort("{normal_cn} column in {.arg mut_data} or {.arg cnv_data} must be numeric")
            }
            data.table::setnames(out, normal_cn, "normal_cn")
        }
    }
    columns <- c(
        on_patient, on_sample, "purity", start_field, end_field,
        on_chr, mut_pos, "alt_counts", "ref_counts"
    )
    if (any(ccf_type == c("phyloCCF", "both"))) {
        data.table::setnames(
            out, c("nAraw", "nBraw"),
            c("major_raw", "minor_raw")
        )
        out <- estimate_phylo_ccf(out,
            patient_field = on_patient,
            sample_field = on_sample,
            chr_field = on_chr,
            pos_field = mut_pos,
            ref_field = ref_field,
            alt_field = alt_field, indel_field = indel_field,
            min_vaf_to_explain = min_vaf_to_explain,
            min_subclonal = min_subclonal, conipher = conipher
        )
        columns <- c(columns, "minor_cn", "major_cn", ref_field, alt_field)
    }
    if (any(ccf_type == c("bootstrapCCF", "both"))) {
        assert_pkg("sequenza")
        assert_pkg("boot")
        data.table::setnames(
            out, c("nMinor", "nMajor"),
            c("minor_cn", "major_cn")
        )
        # here will modify data in place
        estimate_btstr_ccf(out, sample_field = group)
        columns <- c(columns, "minor_cn", "major_cn")
    }
    columns <- c(
        columns, "normal_cn", "expVAF", "obsVAF",
        "mut_multi", "mut_multi_lower", "mut_multi_higher",
        "mut_multi_btstr_lower", "mut_multi_btstr_higher",
        "CCF", "CCF_lower", "CCF_higher",
        "CCF_btstr_lower", "CCF_btstr_higher",
        "absCCF", "absCCF_lower", "absCCF_higher", "phyloCCF",
        "phyloCCF_lower", "phyloCCF_higher", "mutCopyNum",
        "no.chrs.bearing.mut", "whichFrac", "CPNChange",
        kept_cols
    )
    out[, .SD, .SDcols = intersect(names(out), columns)]
}

define_normal_cn <- function(gender, chr) {
    allosomes <- GenomeInfoDb::seqlevelsInGroup(chr, group = "sex")
    autosomes <- GenomeInfoDb::seqlevelsInGroup(chr, group = "auto")
    if (!all(as.character(chr) %chin% c(allosomes, autosomes))) {
        cli::cli_abort(c(
            "Find chromosomes not in {.filed allosomes} and {.field autosomes}",
            i = "Please check {.code ?GenomeInfoDb::seqlevelsInGroup} for definition of {.filed allosomes} and {.field autosomes}"
        ))
    }
    data.table::fifelse(
        gender == "male" & as.character(chr) %chin% allosomes,
        1L, 2L
    )
}
