#' Matching mutation Copy number and Estimating CCF
#'
#' @description
#' Following methods to define subclone and clone.
#' - define clone and subclone based on absCCF or mut_multi_btstr
#' - define early or late with Mt or phyloCCF.
#'
#' For CONIPHER anlayis, use conipher = TRUE, min_vaf_to_explain = 0.05. And you
#' should provide indel_field or both ref_field and alt_field.
#'
#' @param mut_data A data.frame containing mutation data. In addition to the
#' columns defined in the arguments on_patient, on_sample, on_chr, and mut_pos,
#' it is also required that the columns "ref_counts" and "alt_counts" are
#' included.
#' @inheritParams identify_mut_cn
#' @param ccf_type which CCF should we estimate, one of "phyloCCF" (include both
#' absolute and phylogenetic CCF) or "bootstrapCCF". You can provide
#' c("phyloCCF", "bootstrapCCF") to calculate both.
#' @param subclonal_cn_correction A logical indicates whether correcting
#'  subclonal copy number. Only used when phyloCCF is calculated. If TRUE, the
#'  raw copy number (point number) ("nAraw" and "nBraw") will be required.
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
#'   NULL, all chromosome will be used.
#' @param ... Not used currently
#' @param ref_field,alt_field A string specifying the column of reference allele
#'  or variant allele in mut_data to estimate whether a variant is a indel or
#'  not. `alt_field` will only be used when `indel_field` is NULL. `ref_field`
#'  must be specified when `indel_field` or `alt_field` is not NULL.
#' @param indel_field A string specifying a logical column of mut_data which
#' indicates whether or not a variant is a indel.
#' @param min_vaf_to_explain A numeric, the minimal vaf value to define
#'  subclone.
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
    ccf_type = "phyloCCF", subclonal_cn_correction = TRUE,
    on_patient = NULL, on_sample = NULL,
    purity_field = "purity", on_chr = "chr", mut_pos = "pos",
    start_field = "startpos", end_field = "endpos",
    normal_cn = NULL, gender_field = "gender", contigs = NULL,
    ...,
    # arguments for ccf_type = phyloCCF
    ref_field = NULL,
    alt_field = NULL, indel_field = NULL,
    min_vaf_to_explain = NULL, conipher = FALSE,
    # Other arguments
    nomatch = NULL, kept_cols = NULL) {
    assert_pkg("GenomeInfoDb")
    rlang::check_dots_empty()

    # check arguments firstly
    assert_in(ccf_type, c("phyloCCF", "bootstrapCCF"))
    assert_(subclonal_cn_correction, rlang::is_scalar_logical,
        "a scalar {.cls logical}",
        null_ok = TRUE, show_length = TRUE
    )
    assert_(on_patient, rlang::is_scalar_character,
        "a scalar {.cls character}",
        null_ok = TRUE, show_length = TRUE
    )
    assert_(on_sample, rlang::is_scalar_character,
        "a scalar {.cls character}",
        null_ok = TRUE, show_length = TRUE
    )
    assert_(purity_field, rlang::is_scalar_character,
        "a scalar {.cls character}",
        null_ok = FALSE, show_length = TRUE
    )
    assert_(
        on_chr, rlang::is_scalar_character,
        "a scalar {.cls character}", show_length = TRUE
    )
    assert_(
        mut_pos, rlang::is_scalar_character,
        "a scalar {.cls character}", show_length = TRUE
    )
    assert_(
        start_field, rlang::is_scalar_character,
        "a scalar {.cls character}", show_length = TRUE
    )
    assert_(
        end_field, rlang::is_scalar_character,
        "a scalar {.cls character}", show_length = TRUE
    )
    assert_(normal_cn,
        function(x) {
            is_scalar_numeric(x) || rlang::is_scalar_character(x)
        }, "a number or a string",
        null_ok = TRUE, show_length = TRUE
    )
    assert_data_frame(mut_data)
    assert_data_frame_columns(mut_data, c(
        names(on_patient) %||% on_patient,
        names(on_sample) %||% on_sample,
        names(on_chr) %||% on_chr, mut_pos,
        "ref_counts", "alt_counts"
    ))
    cnv_columns <- NULL
    if (any(ccf_type == "phyloCCF") && subclonal_cn_correction) {
        cnv_columns <- c(cnv_columns, c("nAraw", "nBraw"))
    }
    if (any(ccf_type == "bootstrapCCF") ||
        (any(ccf_type == "phyloCCF") && !subclonal_cn_correction)) {
        cnv_columns <- c(cnv_columns, c("nMinor", "nMajor"))
    }
    assert_data_frame(cnv_data)
    assert_data_frame_columns(cnv_data, c(
        on_patient, on_sample, on_chr, start_field, end_field, cnv_columns
    ))
    assert_(kept_cols, is.character, "an atomic {.cls character}",
        null_ok = TRUE, show_length = TRUE
    )
    mut_data <- data.table::as.data.table(mut_data)
    cnv_data <- data.table::as.data.table(cnv_data)

    if (any(ccf_type == "phyloCCF")) {
        # for phyloCCF, we need raw copy number (point number)
        if (subclonal_cn_correction) {
            # define subclone copy number
            cnv_data <- define_subclone_cn(cnv_data, min_subclonal = 0.01)
            data.table::setnames(
                cnv_data, c("nAraw", "nBraw"),
                c("major_raw", "minor_raw")
            )
        } else {
            cnv_data[, c("fracA", "fracB", "fracC", "fracD", "nMaj_A", "nMaj_B", "nMaj_C", "nMaj_D", "nMin_A", "nMin_B", "nMin_C", "nMin_D") := list( # nolint
                1L, 0L, NA_real_, NA_real_,
                nMajor, nMajor, NA_integer_, NA_integer_,
                nMinor, nMinor, NA_integer_, NA_integer_
            )]
            cnv_data[, c("minor_raw", "major_raw") := list(
                nMinor, nMajor
            )]
        }
    }
    if (any(ccf_type == "bootstrapCCF")) {
        data.table::setnames(
            cnv_data, c("nMinor", "nMajor"),
            c("minor_cn", "major_cn")
        )
    }

    # just extract the segmented CNV for this sample
    out <- mut_match_cn(mut_data, cnv_data,
        on_patient = on_patient, on_sample = on_sample, on_chr = on_chr,
        mut_pos = mut_pos, start_field = start_field,
        end_field = end_field, nomatch = nomatch
    )

    # assert every sample has only one patient
    if (!is.null(on_sample) && !is.null(on_patient)) {
        assert_nest(out, on_patient, on_sample)
        init_msg <- "Processing {.val {nrow(out)}} mutation{?s} of {.val {length(unique(out[[on_sample]]))}} sample{?s} from {.val {length(unique(out[[on_patient]]))}} patient{?s}"
    } else if (!is.null(on_patient)) {
        group <- on_patient
        init_msg <- "Processing {.val {nrow(out)}} mutation{?s} from {.val {length(unique(out[[on_patient]]))}} patient{?s} (each has one sample)"
    } else if (!is.null(on_sample)) {
        group <- on_sample
        init_msg <- "Processing {.val {nrow(out)}} mutation{?s} of {.val {length(unique(out[[on_sample]]))}} sample{?s} from {.val {1}} patient"
    } else {
        group <- NULL
        init_msg <- "Processing {.val 1} sample"
    }

    if (is.null(group)) {
        info_msg <- "try to set {.arg on_patient} or {.arg on_sample}"
    } else {
        info_msg <- NULL
    }

    # assert every necessary column is in the combined data
    assert_data_frame_columns(out, c(purity_field, kept_cols),
        arg = c("mut_data", "cnv_data")
    )

    # assert every samples provided only have one purity value
    assert_nest(out, purity_field, group, info_msg = info_msg)
    data.table::setnames(out, purity_field, "purity")
    if (!all(data.table::between(out$purity, 0L, 1L))) {
        cli::cli_abort("purity must in [0, 1]")
    }
    all_seqs <- as.character(out[[on_chr]])
    if (!is.null(contigs)) {
        # filter contigs
        contigs <- as.character(contigs)
        contigs <- map_seqnames(contigs, GenomeInfoDb::seqlevelsStyle(all_seqs))
        matched_contigs <- all_seqs %chin% contigs
        out <- out[matched_contigs]
        all_seqs <- all_seqs[matched_contigs]
    }

    if (nrow(out) == 0L) {
        contig_removing_msg <- "No mutation to proceed"
        if (!is.null(contigs) && !any(matched_contigs)) {
            contig_removing_msg <- paste(
                contig_removing_msg,
                "after filtering by {.arg contigs}"
            )
            contig_removing_msg <- c(contig_removing_msg,
                i = "Please check {.arg contigs}"
            )
        }
        cli::cli_abort(contig_removing_msg)
    }

    allosomes <- GenomeInfoDb::seqlevelsInGroup(all_seqs, group = "sex")
    autosomes <- GenomeInfoDb::seqlevelsInGroup(all_seqs, group = "auto")
    if (!all(all_seqs %in% c(allosomes, autosomes))) {
        cli::cli_abort(c(
            "Find chromosomes not in {.filed allosomes} and {.field autosomes}",
            i = "Please check {.code ?GenomeInfoDb::seqlevelsInGroup} for definition of {.filed allosomes} and {.field autosomes}"
        ))
    }

    # define normal_cn
    if (is.null(normal_cn)) {
        cli::cli_inform(
            "Using {.arg gender_field} to define {.field normal_cn}"
        )
        assert_(
            gender_field, rlang::is_scalar_character,
            "a scalar {.cls character}", show_length = TRUE
        )
        # assert every samples provided only one gender value
        assert_data_frame_columns(out, gender_field,
            arg = c("mut_data", "cnv_data")
        )
        if (!all(out[[gender_field]] %in% c("male", "female"))) {
            cli::cli_abort("Only {.val male} and {.val female} are supported in {.field {gender_field}} column")
        }
        assert_nest(
            out, gender_field, group,
            cross_format = "group", info_msg = info_msg
        )
        out$normal_cn <- define_normal_cn(
            out[[gender_field]], all_seqs, allosomes
        )
    } else {
        if (is_scalar_numeric(normal_cn)) {
            out[, normal_cn := normal_cn]
        } else {
            assert_data_frame_columns(out, normal_cn,
                arg = "mut_data"
            )
            if (!is.numeric(out[[normal_cn]])) {
                cli::cli_abort("{normal_cn} column in {.arg mut_data} or {.arg cnv_data} must be numeric")
            }
            data.table::setnames(out, normal_cn, "normal_cn")
        }
    }
    cli::cli_inform(init_msg)
    columns <- c(
        on_patient, on_sample, "purity", start_field, end_field,
        on_chr, mut_pos, "alt_counts", "ref_counts"
    )
    if (any(ccf_type == "phyloCCF")) {
        out <- estimate_phylo_ccf(out,
            patient_field = on_patient,
            sample_field = on_sample,
            chr_field = on_chr,
            pos_field = mut_pos,
            ref_field = ref_field,
            alt_field = alt_field, indel_field = indel_field,
            min_vaf_to_explain = min_vaf_to_explain,
            conipher = conipher
        )
        columns <- c(columns, ref_field, alt_field)
    }
    if (any(ccf_type == "bootstrapCCF")) {
        # here will modify data in place
        estimate_btstr_ccf(out, sample_field = group)
    }
    columns <- c(
        columns, "minor_raw", "major_raw", "minor_cn", "major_cn",
        "normal_cn", "expVAF", "obsVAF",
        "mut_multi", "mut_multi_lower", "mut_multi_higher",
        "mut_multi_btstr_lower", "mut_multi_btstr_higher",
        "CCF", "CCF_lower", "CCF_higher",
        "CCF_btstr_lower", "CCF_btstr_higher",
        "absCCF", "absCCF_lower", "absCCF_higher", "phyloCCF",
        "phyloCCF_lower", "phyloCCF_higher", "mutCopyNum",
        "no.chrs.bearing.mut", "whichFrac", "CPNChange",
        kept_cols
    )
    out[, .SD, .SDcols = intersect(columns, names(out))]
}

define_normal_cn <- function(gender, chr, allosomes = NULL) {
    allosomes <- allosomes %||%
        GenomeInfoDb::seqlevelsInGroup(chr, group = "sex")
    data.table::fifelse(gender == "male" & chr %in% allosomes, 1L, 2L)
}
