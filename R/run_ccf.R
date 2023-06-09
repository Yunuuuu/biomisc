#' Matching mutation Copy number and Estimating CCF (absolute and phylogenetic)
#' @description
#' - define clone and subclone based on absCCF or mut.multi.bstr
#' - define early or late with Mt or phyloCCF.
#'
#' For CONIPHER anlayis, use min_subclonal = 0.05, conipher = TRUE,
#' min_vaf_to_explain = 0.05. And you should provide indel_field or both
#' ref_field and alt_field.
#'
#' @inheritParams identify_mut_cn
#' @param purity_field A string specifying the purity column (in `mut_data` or
#'  `cnv_data`). Default is "purity".
#' @param indel_field A string specifying a logical column of mut_cn_data which
#' indicates whether or not a variant is a indel for indel phyloCCF corretion.
#' @param ref_field,alt_field A string specifying the column of reference allele
#'  or variant allele in mut_cn_data to estimate whether a variant is a indel or
#'  not. This will only be used when indel_field is NULL.
#' @inheritDotParams estimate_ccf -mut_cn_data -sample_field -chr_field
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
#' @inherit estimate_ccf note
#' @export
run_ccf <- function(
    mut_data, cnv_data, on_sample = NULL, purity_field = NULL,
    on_chr = "chr", mut_pos = "pos", start_field = "startpos",
    end_field = "endpos", indel_field = NULL, ref_field = NULL,
    alt_field = NULL, ..., nomatch = NULL, kept_cols = NULL) {
    assert_class(on_sample, rlang::is_scalar_character,
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
    assert_df_with_columns(mut_data, c(
        names(on_sample) %||% on_sample,
        names(on_chr) %||% on_chr, mut_pos,
        "ref_counts", "alt_counts"
    ))
    assert_df_with_columns(cnv_data, c(
        on_sample, on_chr, start_field, end_field,
        "nAraw", "nBraw"
    ))
    assert_class(kept_cols, is.character, "{.cls character} vector",
        cross_msg = NULL, null_ok = TRUE
    )

    mut_data <- data.table::as.data.table(mut_data)
    cnv_data <- data.table::as.data.table(cnv_data)

    # define subclone copy number
    cnv_data <- define_subclone_cn(cnv_data, min_subclonal = 0.01)

    # just extract the segmented CNV for this sample
    out <- mut_match_cn(mut_data, cnv_data,
        on_sample = on_sample, on_chr = on_chr,
        mut_pos = mut_pos, start_field = start_field,
        end_field = end_field, nomatch = nomatch
    )
    purity_field <- purity_field %||% "purity"
    assert_df_with_columns(out, c(purity_field, kept_cols),
        check_class = FALSE,
        arg = c("mut_data", "cnv_data")
    )
    data.table::setnames(out, c("nAraw", "nBraw"), c("major_raw", "minor_raw"))
    out <- estimate_ccf(out,
        sample_field = on_sample,
        purity_field = purity_field,
        chr_field = on_chr, ...
    )
    # to ensure potentially unreliable VAFs of indels did not lead to
    # separate mutation clusters, each estimated indel CCF was multiplied by a
    # region specific correction factor.  Assuming the majority of ubiquitous
    # mutations, present in all regions, are clonal, the region specific
    # correction factor was calculated by dividing the median mutation CCF of
    # ubiquitous mutations by the median indel CCF of ubiquitous indels.

    # SNV: width(ref) == width(alt) including snv, dbs, mbs
    if (is.null(indel_field)) {
        if (!is.null(ref_field) && !is.null(alt_field)) {
            out$..tmp_is_indel.. <- is_indel(
                out[[ref_field]], out[[alt_field]]
            )
        } else if (!is.null(ref_field) || !is.null(alt_field)) {
            cli::cli_abort(
                "both {.arg ref_field} and {.arg alt_field} must be specified to correct indel phyloCCF"
            )
        }
    } else {
        if (!is.logical(out[[indel_field]])) {
            cli::cli_abort("{indel_field} column in {.arg out} must be logical")
        }
        data.table::setnames(out, indel_field, "..tmp_is_indel..")
    }

    if (!is.null(out$..tmp_is_indel..)) {
        if (is.null(on_sample)) {
            ..tmp_nsamples.. <- 1L
        } else {
            ..tmp_nsamples.. <- length(unique(out[[on_sample]]))
        }
        # nolint start
        ubiq_out <- out[, .SD[sum(obsVAF > 5L) == ..tmp_nsamples..],
            by = c(on_chr, mut_pos, ref_field)
        ][!is.na(phyloCCF)]
        ..tmp_indel_cf.. <- ubiq_out[, list(
            ..tmp_correction_factor =
                median(phyloCCF[!..tmp_is_indel..]) /
                    median(phyloCCF[..tmp_is_indel..])
        ), by = on_sample]
        if (is.null(on_sample)) {
            ..tmp_indel_cf.. <- ..tmp_indel_cf..$..tmp_correction_factor
        } else {
            ..tmp_indel_cf.. <- structure(
                ..tmp_indel_cf..$..tmp_correction_factor,
                names = ..tmp_indel_cf..[[on_sample]]
            )
            ..tmp_indel_cf.. <- ..tmp_indel_cf..[out[[on_sample]]]
        }
        out[(..tmp_is_indel..), phyloCCF := phyloCCF * ..tmp_indel_cf..]
        out[(..tmp_is_indel..), mutCopyNum := mutCopyNum * ..tmp_indel_cf..]
        out[, ..tmp_is_indel.. := NULL]
        # nolint end
    }

    # out[, c("major_cn", "minor_cn") := list(
    #     major_cn = pmax(nMinor, nMajor), # nolint
    #     minor_cn = pmin(nMinor, nMajor) # nolint
    # )]
    columns <- c(
        on_sample, purity_field, on_chr, mut_pos, ref_field, alt_field,
        "alt_counts", "ref_counts", start_field, end_field
    )
    for (i in c("nMajor", "nMinor")) {
        if (any(i == names(cnv_data))) {
            columns <- c(columns, i)
        }
    }
    columns <- c(
        columns, "major_raw", "minor_raw",
        "nMaj_A", "nMaj_B", "nMaj_C", "nMaj_D",
        "nMin_A", "nMin_B", "nMin_C", "nMin_D",
        "fracA", "fracB", "fracC", "fracD",
        "normal_cn", "expVAF", "obsVAF", "absCCF",
        "absCCF_lower", "absCCF_higher", "phyloCCF",
        "phyloCCF_lower", "phyloCCF_higher", "mutCopyNum",
        "no.chrs.bearing.mut", "whichFrac", "CPNChange",
        kept_cols
    )
    out[, .SD, .SDcols = columns]
}

#' Estimating CCF (absolute and phylogenetic)
#'
#' @param mut_cn_data A data.frame with mutation and copy number data. Copy
#'  number often contain subclonal copy number as described in
#'  [CONIPHER](https://github.com/McGranahanLab/CONIPHER-wrapper/blob/b58235d1cb42d5c7fd54122dc6b9f5e6c4110a75/src/TRACERxHelperFunctions.R#L1).
#' @param sample_field A string specifying the sample column. If NULL, all data
#'  will be regarded from the same sample. This is used to confirm every sample
#'  have the same `purity`.
#' @param purity_field A string specifying the purity column. Default is
#'  "purity".
#' @param contigs An atomic vector specifying the chromosome to analyze.
#' @param chr_field A string specifying the chromosome column. Only used when
#'   `contigs` is not NULL or `normal_cn` is NULL. Default is "chr".
#' @param normal_cn A scalar number specifying the normal.copy number or a
#'  string indicating the normal_cn column in `mut_cn_data`.  It's save to use 2
#'  if you only analyze autosomes. Or you should use `gender` to define the
#'  `normal_cn`. For sex chromosomes and gender is "male", 1L will be used,
#'  otherwise, 2L will be used.
#' @param gender A string, only "female" and "male" are allowed.
#' @param min_subclonal Minimal copy number to define subclone.
#' @param conipher A scalar logical indicates whether calculate phyloCCF like
#'  CONIPHER. Details see
#'  <https://github.com/McGranahanLab/CONIPHER-wrapper/tree/main>
#' @param min_vaf_to_explain A numeric, the minimal vaf value to define
#'  subclone.
#' @note This function should run for every patient of all region samples.
#' @seealso [run_ccf]
#' @export
estimate_ccf <- function(mut_cn_data, sample_field = NULL, purity_field = NULL, contigs = NULL, chr_field = NULL, normal_cn = 2L, gender, min_vaf_to_explain = NULL, min_subclonal = NULL, conipher = FALSE) {
    # check arguments firstly
    assert_class(purity_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE,
        cross_msg = NULL
    )
    purity_field <- purity_field %||% "purity"
    assert_df_with_columns(mut_cn_data, c(
        sample_field, purity_field,
        "major_raw", "minor_raw", "alt_counts", "ref_counts",
        "nMaj_A", "nMaj_B", "nMaj_C", "nMaj_D",
        "nMin_A", "nMin_B", "nMin_C", "nMin_D",
        "fracA", "fracB", "fracC", "fracD"
    ))
    assert_class(sample_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE,
        cross_msg = NULL
    )
    assert_class(chr_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE,
        cross_msg = NULL
    )
    assert_class(normal_cn,
        function(x) {
            is_scalar_numeric(x) || rlang::is_scalar_character(x)
        }, "scalar {.cls numeric} or {.cls character}",
        null_ok = TRUE, cross_msg = NULL
    )
    assert_class(min_subclonal, is_scalar_numeric,
        "scalar {.cls numeric}",
        null_ok = TRUE, cross_msg = NULL
    )
    assert_class(conipher, rlang::is_scalar_logical,
        "scalar {.cls logical}",
        null_ok = TRUE, cross_msg = NULL
    )
    mut_cn_data <- data.table::as.data.table(mut_cn_data)
    if (!is.null(contigs) || is.null(normal_cn)) {
        chr_field <- chr_field %||% "chr"
        assert_df_with_columns(mut_cn_data, chr_field, check_class = FALSE)
    }
    if (!is.null(contigs)) {
        # filter contigs
        matched_contigs <- as.character(mut_cn_data[[chr_field]]) %chin%
            as.character(contigs)
        mut_cn_data <- mut_cn_data[matched_contigs]
    }
    if (is.null(normal_cn)) {
        # assert every patient provided only one gender value
        gender <- match.arg(tolower(gender), c("female", "male"))
        mut_cn_data$normal_cn <- define_normal_cn(
            gender, mut_cn_data[[chr_field]]
        )
    } else {
        if (is_scalar_numeric(normal_cn)) {
            mut_cn_data$normal_cn <- normal_cn
        } else {
            assert_df_with_columns(mut_cn_data, normal_cn, check_class = FALSE)
            if (!is.numeric(mut_cn_data[[normal_cn]])) {
                cli::cli_abort("{normal_cn} column in {.arg mut_cn_data} must be numeric")
            }
            data.table::setnames(mut_cn_data, normal_cn, "normal_cn")
        }
    }
    # assert every samples provided only have one purity value
    if (!all(data.table::between(mut_cn_data[[purity_field]], 0L, 1L))) {
        cli::cli_abort("purity must in [0, 1]")
    }
    assert_nest(
        mut_cn_data, purity_field, sample_field,
        info_msg = if (is.null(sample_field)) "try to set {.arg sample_field}"
    )

    # In order to estimate whether mutations were clonal or subclonal, and the
    # clonal structure of each tumor, a modified version of PyClone was used.
    # For each mutation, two values were calculated, obsCCF and phyloCCF. obsCCF
    # corresponds to the observed cancer cell fraction (CCF) of each mutation.
    # Conversely, phyloCCF corresponds to the phylogenetic CCF of a mutation. To
    # clarify the difference between these two values, consider a mutation
    # present in every cancer cell within a tumor. A subclonal copy number event
    # in one tumor region may lead to loss of this mutation in a subset of
    # cancer cells.  While, the obsCCF of this mutation is therefore below 1,
    # from a phylogenetic perspective the mutation can be considered ‘clonal’ as
    # it occurred on the trunk of the tumor’s phylogenetic tree, and, as such,
    # the phyloCCF may be 1.
    out <- mut_cn_data

    # define normal_cn
    # nolint start
    out[, expVAF := calculate_vaf(
        1L, purity, major_raw + minor_raw, normal_cn
    )]
    out[, obsVAF := alt_counts / (alt_counts + ref_counts)]
    absolute_ccfs <- calculate_abs_ccf(
        out$alt_counts, out$ref_counts,
        CNts = out$major_raw + out$minor_raw,
        purity = out$purity, CNns = out$normal_cn
    )
    out$absCCF <- absolute_ccfs$est
    out$absCCF_lower <- absolute_ccfs$lower
    out$absCCF_higher <- absolute_ccfs$higher

    out[
        ,
        c("phyloCCF", "phyloCCF_lower", "phyloCCF_higher") := calculate_ccf(
            alt_counts, ref_counts,
            CNts = major_raw + minor_raw,
            purity, observed_vafs = obsVAF,
            expected_vafs = expVAF, CNns = normal_cn
        )
    ]
    # nolint end
    out$mutCopyNum <- out$phyloCCF
    # pre-allocate column for downstream analysis
    out$no.chrs.bearing.mut <- 1
    out$whichFrac <- NA_character_
    out$CPNChange <- 0L

    # Identification of subclonal mutations  -----------------------
    # nolint start
    out[, ..tmp_is_subclone.. := mutCopyNum > 0.01]
    if (!is.null(min_vaf_to_explain)) {
        out[, ..tmp_is_subclone.. := ..tmp_is_subclone.. & obsVAF >= min_vaf_to_explain]
    }
    if (conipher) {
        # if the observed variant allele frequency was significantly different
        # from that expected (P<0.01, using prop.test in R) given a clonal
        # mutation, we determined whether a subclonal copy number event could
        # result in a non-significant (P>0.01) difference between observed and
        # expected VAFs.
        out[, ..tmp_is_subclone.. := ..tmp_is_subclone.. & prop_test_pvalues(
            alt_counts, alt_counts + ref_counts, trim_value(expVAF),
            alternative = "less",
            correction = "BH"
        ) < 0.01]
        out[, ..tmp_is_subclone.. := ..tmp_is_subclone.. & absCCF_higher < 1L]
    } else {
        out[
            ,
            ..tmp_is_subclone.. := ..tmp_is_subclone.. & absolute_ccfs$prob.subclonal > 0.5
        ]
    }

    out[..tmp_is_subclone.. & fracA == 1L, whichFrac := "A,B"]
    # nolint end

    # define best CN -----------------------------------------------
    # nolint start
    out[
        ..tmp_is_subclone.. & fracA != 1L & is.na(fracC),
        ..tmp_best_cn.. := calculate_best_cn_for_loss_mut(
            nMaj_A, nMaj_B, nMin_A, nMin_B,
            fracA, fracB, fracA, fracB, mutCopyNum
        )
    ]
    out[
        ..tmp_is_subclone.. & fracA != 1L & !is.na(fracC),
        ..tmp_best_cn.. := calculate_best_cn_for_loss_mut(
            nMaj_A + nMaj_B, nMaj_C + nMaj_D,
            nMin_A + nMin_D, nMin_C + nMin_B,
            fracA + fracB, fracC + fracD,
            fracA + fracD, fracC + fracB, mutCopyNum
        )
    ]
    out[, ..tmp_matched_rows.. := ..tmp_best_cn.. == 1L]
    if (!is.null(min_subclonal)) {
        out[
            ,
            ..tmp_matched_rows.. := ..tmp_matched_rows.. |
                ..tmp_best_cn.. >= (1 - min_subclonal)
        ]
    }
    out[
        (..tmp_matched_rows..),
        whichFrac := data.table::fifelse( # nolint
            is.na(fracC), "A,B", "A,B,C,D" # nolint
        )
    ]
    # nolint end

    # check whether subclonal CN results in clonal mutation
    # otherwise subclonal CN doesn't explain subclonal MCN
    # nolint start
    out[, ..tmp_expProp.. := expVAF * ..tmp_best_cn..]
    out[
        !(..tmp_matched_rows..), # nolint
        ..tmp_explained_by_cn_pvalue.. := prop_test_pvalues( # nolint
            alt_counts, (alt_counts + ref_counts) * purity, # nolint
            prop = ..tmp_expProp.. / purity, # nolint
            alternative = "less"
        )
    ]
    out[, ..tmp_matched_rows.. := NULL] # nolint
    if (conipher) {
        out[
            ..tmp_explained_by_cn_pvalue.. > 0.01,
            ..tmp_operated_rows.. := calculate_mut_multi(
                major_raw + minor_raw,
                prop_test_ci(
                    alt_counts, alt_counts + ref_counts,
                    trim_value(..tmp_expProp..)
                )[[1L]],
                purity
            ) / ..tmp_best_cn.. <= 1L
        ]
    } else {
        out[, ..tmp_operated_rows.. := ..tmp_explained_by_cn_pvalue.. > 0.01]
    }
    # nolint end
    out[
        (..tmp_operated_rows..), # nolint
        c("phyloCCF", "phyloCCF_lower", "phyloCCF_higher", "no.chrs.bearing.mut", "expVAF", "CPNChange") := {
            tmp_vafs <- prop_test_ci(
                alt_counts, alt_counts + ref_counts, # nolint
                trim_value(..tmp_expProp..) # nolint
            )
            tmp_CNts <- major_raw + minor_raw # nolint
            list(
                mutCopyNum / ..tmp_best_cn.., # nolint
                calculate_mut_multi(
                    CNts = tmp_CNts,
                    vafs = tmp_vafs[[1L]],
                    purity = purity # nolint
                ) / ..tmp_best_cn..,
                calculate_mut_multi(
                    CNts = tmp_CNts,
                    vafs = tmp_vafs[[2L]],
                    purity = purity,
                ) / ..tmp_best_cn..,
                ..tmp_best_cn.., ..tmp_expProp.., 1L # nolint
            )
        }
    ]
    out[, ..tmp_explained_by_cn_pvalue.. := NULL] # nolint
    out[, ..tmp_mut_multi.. := NULL] # nolint
    out[, ..tmp_expProp.. := NULL] # nolint
    out[, ..tmp_operated_rows.. := NULL] # nolint
    out[, ..tmp_is_subclone.. := NULL] # nolint
    out[, ..tmp_best_cn.. := NULL] # nolint

    # Next, let's deal with potentially amplified mutations
    # convert MCN to subclonal fraction - tricky for amplified mutations test
    # for mutations in more than 1 copy
    out[
        ,
        ..tmp_amp_mut_pvalue.. := prop_test_pvalues( # nolint
            alt_counts, alt_counts + ref_counts, trim_value(expVAF), # nolint
            alternative = "greater"
        )
    ]

    # copy numbers of subclones can only differ by 1 or 0 (as assumed when
    # calling subclones)
    # nolint start
    out[..tmp_amp_mut_pvalue.. <= 0.05 & mutCopyNum > 1L & is.na(fracC), c("no.chrs.bearing.mut", "..tmp_best_cn..") := {
        calculate_best_cn_for_amp_mut(
            nMaj_A, nMaj_B,
            fracA = fracA, fracB = fracB,
            alt_counts = alt_counts, mutCopyNum = mutCopyNum,
            conipher = conipher
        )
    }]
    out[..tmp_amp_mut_pvalue.. <= 0.05 & mutCopyNum > 1L & !is.na(fracC), c("no.chrs.bearing.mut", "..tmp_best_cn..") := {
        calculate_best_cn_for_amp_mut(
            nMaj_A, nMaj_C, nMin_A, nMin_B,
            fracA + fracB, fracC + fracD,
            fracA + fracD, fracC + fracB,
            alt_counts = alt_counts, mutCopyNum = mutCopyNum,
            conipher = conipher
        )
    }]
    out[..tmp_amp_mut_pvalue.. <= 0.05 & mutCopyNum > 1L, c("phyloCCF", "phyloCCF_lower", "phyloCCF_higher", "expVAF") := {
        list(
            mutCopyNum / ..tmp_best_cn..,
            phyloCCF_lower / ..tmp_best_cn..,
            phyloCCF_higher / ..tmp_best_cn..,
            expVAF * ..tmp_best_cn..
        )
    }]
    out[, ..tmp_amp_mut_pvalue.. := NULL]
    out[, ..tmp_best_cn.. := NULL]
    # nolint end

    # finally, let's sort out 'missing' ones
    out[
        alt_counts == 0L, # nolint
        c("no.chrs.bearing.mut", "expVAF", "absCCF") := list(0L, 0L, 0L)
    ]
    out[]
}

#' min.subclonal was set to 0.1 in
#' https://bitbucket.org/nmcgranahan/clonalneoantigenanalysispipeline/src/master/
#' @return A data.table
#' @keywords validated
#' @noRd
define_subclone_cn <- function(seg, min_subclonal = 0.01) {
    seg <- data.table::as.data.table(seg)
    seg[, c("nMaj1", "nMin1") := lapply(.SD, floor),
        .SDcols = c("nAraw", "nBraw")
    ]
    seg[, c("nMaj2", "nMin2") := lapply(.SD, ceiling),
        .SDcols = c("nAraw", "nBraw")
    ]
    seg[, fracMaj1 := nMaj2 - nAraw] # nolint
    seg[, fracMin1 := nMin2 - nBraw] # nolint
    seg[, c("fracMaj2", "fracMin2") := lapply(.SD, function(x) {
        1L - x
    }), .SDcols = c("fracMaj1", "fracMin1")]

    # how much of the genome for each tumour region is subject to subclonal copy
    # number
    # prop.aber <- seg[, # nolint
    #     {
    #         prop_maj1 <- fracMaj1 < min_subclonal | # nolint
    #             fracMaj1 > 1 - min_subclonal
    #         prop_min1 <- fracMin1 < min_subclonal | # nolint
    #             fracMin1 > 1 - min_subclonal
    #         sum(prop_maj1 & prop_min1, na.rm = TRUE) / .N
    #     },
    #     by = c(names(seg)[[1L]])
    # ] # the original function return a atomic vector

    # next, let's deal with the minimum subclonal
    seg[, c("fracMaj1", "fracMaj2", "fracMin1", "fracMin2") := lapply(.SD, function(x) {
        data.table::fcase(
            x < min_subclonal, 0,
            x > 1 - min_subclonal, 1,
            rep_len(TRUE, length(x)), x
        )
    }), .SDcols = c("fracMaj1", "fracMaj2", "fracMin1", "fracMin2")]

    # which ones work?
    sorted_idx <- (seg$fracMaj1 == seg$fracMin1 |
        seg$fracMaj1 == seg$fracMin2) |
        seg$fracMaj1 == 1L |
        seg$fracMaj1 == 0L |
        seg$fracMin1 == 1L |
        seg$fracMin1 == 0L
    seg_sorted <- seg[sorted_idx]
    seg_problem <- seg[!sorted_idx]

    # let's divide the problem further
    seg_problem[, fracA := fracMaj1 * fracMin1] # nolint
    seg_problem[, fracB := fracMaj1 - fracA] # nolint
    seg_problem[, fracC := fracMaj2 * fracMin2] # nolint
    seg_problem[, fracD := fracMaj2 - fracC] # nolint
    seg_problem[, c("nMaj_A", "nMin_A", "nMaj_B", "nMin_B", "nMaj_C", "nMin_C", "nMaj_D", "nMin_D") := list(
        nMaj1, nMin1, nMaj1, nMin2, nMaj2, nMin2, nMaj2, nMin1 # nolint
    )]

    # let's make the sorted easier to read
    # seg_sorted2 <- data.table::copy(seg_sorted)

    seg_sorted[, c("fracA", "fracB", "fracC", "fracD", "nMaj_A", "nMaj_B", "nMaj_C", "nMaj_D", "nMin_A", "nMin_B", "nMin_C", "nMin_D") := {
        fracA.major <- pmax(fracMaj1, fracMaj2, na.rm = TRUE) # nolint
        fracB.major <- pmin(fracMaj1, fracMaj2, na.rm = TRUE) # nolint
        fracA.minor <- pmax(fracMin1, fracMin2, na.rm = TRUE) # nolint
        fracB.minor <- pmin(fracMin1, fracMin2, na.rm = TRUE) # nolint

        nMaj_A <- data.table::fifelse(
            fracA.major == fracMaj1, nMaj1, nMaj2 # nolint
        )
        nMaj_B <- data.table::fifelse(
            fracA.major == fracMaj1, nMaj2, nMaj1 # nolint
        )
        nMin_A <- data.table::fifelse(
            fracA.minor == fracMin1, nMin1, nMin2 # nolint
        )
        nMin_B <- data.table::fifelse(
            fracA.minor == fracMin1, nMin2, nMin1 # nolint
        )
        fracA <- data.table::fcase(
            fracA.minor == 1L, fracA.major,
            fracA.major == 1L, fracA.minor,
            rep_len(TRUE, .N), fracA.major
        )
        fracB <- data.table::fcase(
            fracA.minor == 1L, fracB.major,
            fracA.major == 1L, fracB.minor,
            rep_len(TRUE, .N), fracB.major
        )
        nMaj_B <- data.table::fifelse(fracA.major == 1L, nMaj_A, nMaj_B)
        nMin_B <- data.table::fifelse(fracA.minor == 1L, nMin_A, nMin_B)
        tmp <- rep_len(NA, .N)
        list(
            fracA, fracB, tmp, tmp,
            nMaj_A, nMaj_B, tmp, tmp,
            nMin_A, nMin_B, tmp, tmp
        )
    }]

    data.table::rbindlist(list(seg_problem, seg_sorted),
        use.names = TRUE, fill = TRUE
    )

    # finally, let's choose the columns we want and the order we want
    # columns <- c(
    #     "SampleID", "chr", "startpos", "endpos", "n.het", "cnTotal", "nMajor",
    #     "nMinor", "Ploidy", "ACF", "nAraw", "nBraw", "fracA", "nMaj_A",
    #     "nMin_A", "fracB", "nMaj_B", "nMin_B", "fracC", "nMaj_C", "nMin_C",
    #     "fracD", "nMaj_D", "nMin_D"
    # )
    # let's order this correctly
    # seg_out[order(SampleID, chr, startpos), .SD, .SDcols = columns] # nolint
}

utils::globalVariables(c(
    "sample_id", "..tmp_amp_mut_pvalue..", "best_cn", "expVAF", "obsVAF", "purity",
    "..tmp_explained_by_cn_pvalue..",
    "fracA", "fracB", "fracC", "fracD",
    "fracMaj1", "fracMaj2", "fracMin1", "fracMin2", "..tmp_is_subclone..",
    "major_raw", "minor_raw", "mutCopyNum", "nAraw", "nBraw",
    "nMajor", "nMinor",
    "nMaj1", "nMaj2", "nMaj_A", "nMaj_B",
    "nMaj_C", "nMaj_D", "nMin1", "nMin2",
    "nMin_A", "nMin_B", "nMin_C", "nMin_D",
    "phyloCCF_higher", "phyloCCF_lower", "absCCF_higher", "..tmp_expProp..",
    "ref_counts", "startpos", "alt_counts", "whichFrac",
    "..tmp_operated_rows..", "..tmp_matched_rows..", "normal_cn", "..tmp_mut_multi.."
))

define_normal_cn <- function(gender, chr) {
    allosomes <- GenomeInfoDb::seqlevelsInGroup(chr, group = "sex")
    autosomes <- GenomeInfoDb::seqlevelsInGroup(chr, group = "auto")
    if (!all(as.character(chr) %chin% c(allosomes, autosomes))) {
        cli::cli_abort(c(
            "Find chromosomes not in {.filed allosomes} and {.field autosomes}",
            i = "Please check {.code ?GenomeInfoDb::seqlevelsInGroup} for definition of {.filed allosomes} and {.field autosomes}"
        ))
    }
    if (gender == "male") {
        data.table::fifelse(as.character(chr) %chin% allosomes, 1L, 2L)
    } else {
        2L
    }
}

# Multiplicity of a mutation: the number of DNA copies bearing a mutation m.
calculate_mut_multi <- function(CNts, vafs, purity, CNns = 2L) {
    (vafs / purity) * ((purity * CNts) + CNns * (1L - purity))
}

# purity: The purity is the fraction of cancer cells in the tumor sample.
# local.copy.number: CNt
# threshold was set to 1 - 1e-6 in
# https://github.com/McGranahanLab/CONIPHER-wrapper/blob/b58235d1cb42d5c7fd54122dc6b9f5e6c4110a75/src/run_clustering.R#L530
# Multiplicity (m): The number of DNA copies bearing a mutation m
calculate_vaf <- function(m, purity, CNts, CNns = 2L, threshold = 1L) {
    out <- (purity * m) / (CNns * (1L - purity) + purity * CNts)
    trim_value(out, threshold = threshold)
}

calculate_abs_ccf <- function(
    alt_counts, ref_counts, CNts, purity,
    CNns = 2L, candidate_mut_multi = seq(0.01, 1, length.out = 100L),
    alpha = 0.05) {
    # use maximum likelihood method to define absCCF
    assert_length(CNns, length = length(alt_counts), scalar_ok = TRUE)
    assert_length(purity, length = length(alt_counts), scalar_ok = TRUE)
    out_list <- .mapply(
        function(alt_count, ref_count, CNt, CNn, p) {
            x <- stats::dbinom(alt_count, alt_count + ref_count,
                prob = calculate_vaf(candidate_mut_multi,
                    purity = p, CNt, CNns = CNn
                )
            )
            if (min(x) == 0L) {
                x[length(x)] <- 1L
            }
            xnorm <- x / sum(x)
            xsort <- sort(xnorm, decreasing = TRUE)
            xcumLik <- cumsum(xsort)
            idx <- xnorm >= xsort[sum(xcumLik < 1L - alpha) + 1L]
            cint <- x[idx]
            mut_multi <- candidate_mut_multi[idx]
            data.table::data.table(
                lower = mut_multi[[1L]],
                est = mut_multi[[which.max(cint)]],
                higher = mut_multi[[length(mut_multi)]],
                prob.subclonal = sum(xnorm[1:90]),
                prob.clonal = sum(xnorm[91:100])
            )
        },
        list(
            alt_count = alt_counts, ref_count = ref_counts,
            CNt = CNts, CNn = rep_len(CNns, length(alt_counts)),
            p = rep_len(purity, length(alt_counts))
        ),
        NULL
    )
    data.table::rbindlist(out_list)
}

calculate_ccf <- function(alt_counts, ref_counts, CNts, purity, observed_vafs = NULL, expected_vafs = NULL, CNns = 2L, prefix = NULL) {
    depths <- alt_counts + ref_counts
    obs_vafs <- observed_vafs %||% (alt_counts / depths)
    if (!is.null(expected_vafs)) {
        expected_vafs <- trim_value(expected_vafs)
    }
    vafs <- prop_test_ci(alt_counts, depths, expected_vafs)
    out_list <- lapply(list(obs_vafs, vafs[[1L]], vafs[[2L]]), function(vaf) {
        calculate_mut_multi(
            CNts = CNts, vafs = vaf,
            purity = purity, CNns = CNns
        )
    })
    if (!is.null(prefix)) {
        names(out_list) <- paste0(prefix, c("", "_lower", "_higher"))
    }
    out_list
}

# check which subclonal mutations can be explained by copy number
prop_test_pvalues <- function(counts, totals, prop, alternative = "two.sided", conf.level = 0.95, correction = NULL) {
    pvalues <- prop_test(counts, totals, prop,
        alternative = alternative,
        conf.level = conf.level
    )[[1L]]
    if (!is.null(correction)) {
        pvalues <- stats::p.adjust(pvalues, method = correction)
    }
    pvalues
}

prop_test_ci <- function(counts, totals, prop, alternative = "two.sided", conf.level = 0.95) {
    prop_test(counts, totals, prop,
        alternative = alternative,
        conf.level = conf.level
    )[2:3]
}

prop_test <- function(counts, totals, prop, alternative = "two.sided", conf.level = 0.95) {
    out_list <- .mapply(
        function(count, total, prop) {
            out <- stats::prop.test(
                count, total, prop,
                alternative = alternative,
                conf.level = conf.level
            )
            c(out$p.value, out$conf.int)
        },
        list(count = counts, total = totals, prop = prop), NULL
    )
    data.table::transpose(out_list)
}

calculate_best_cn_for_loss_mut <- function(
    A, B, C = A, D = B,
    fracA, fracB, fracC = fracA, fracD = fracB,
    mutCopyNum) {
    x <- data.table::fcase(
        A > B, fracA,
        A < B, fracB
    )
    y <- data.table::fcase(
        C > D, fracC,
        C < D, fracD
    )
    possible_frac <- matrix(c(rep_len(1L, length(A)), x, y), ncol = 3L)
    col_idx <- apply(
        abs(mutCopyNum / possible_frac - 1L),
        1L, which.min,
        simplify = TRUE
    )
    possible_frac[cbind(seq_len(nrow(possible_frac)), col_idx)]
}

calculate_best_cn_for_amp_mut <- function(
    A, B, C = NULL, D = NULL,
    fracA, fracB, fracC = NULL, fracD = NULL,
    alt_counts, mutCopyNum, conipher) {
    m_max_cn1 <- pmax(A, B)
    m_frac1_mut <- data.table::fifelse(A < B, fracB, fracA)
    m_max_cn2 <- pmin(A, B)
    m_frac2_mut <- data.table::fifelse(A < B, fracA, fracB)
    arg_list <- list(
        max_cn1 = m_max_cn1, max_cn2 = m_max_cn2,
        frac1_mut = m_frac1_mut, frac2_mut = m_frac2_mut
    )
    if (!is.null(C)) {
        m_max_cn3 <- pmax(C, D)
        m_frac3_mut <- data.table::fifelse(C < D, fracD, fracC)
        m_max_cn4 <- pmin(C, D)
        m_frac4_mut <- data.table::fifelse(C < D, fracC, fracD)
        arg_list <- c(
            arg_list, list(
                max_cn3 = m_max_cn3, max_cn4 = m_max_cn4,
                frac3_mut = m_frac3_mut, frac4_mut = m_frac4_mut
            )
        )
    }
    # I don't know what this does?
    out_list <- .mapply(
        function(max_cn1, max_cn2, frac1_mut, frac2_mut,
                 max_cn3 = NULL, max_cn4 = NULL,
                 frac3_mut = NULL, frac4_mut = NULL,
                 mut_cn, alt_count, conipher = FALSE) {
            best_err <- mut_cn - 1L
            allCNs <- best_cn <- 1L
            for (j in seq_len(max_cn1)) {
                for (k in (j - 1):min(j, max_cn2)) {
                    potential_cn <- j * frac1_mut + k * frac2_mut
                    err <- abs(mut_cn / potential_cn - 1L)
                    if (err < best_err) {
                        best_err <- err
                        best_cn <- potential_cn
                        allCNs <- c(allCNs, best_cn)
                    }
                }
            }
            if (!is.null(max_cn3)) {
                for (j in seq_len(max_cn3)) {
                    for (k in (j - 1):min(j, max_cn4)) {
                        potential_cn <- j * frac3_mut + k * frac4_mut
                        err <- abs(mut_cn / potential_cn - 1L)
                        if (err < best_err) {
                            best_err <- err
                            best_cn <- potential_cn
                            allCNs <- c(allCNs, best_cn)
                        }
                    }
                }
            }
            out <- best_cn # for no.chrs.bearing.mut
            if (conipher) {
                # copied from
                # https://github.com/McGranahanLab/CONIPHER-wrapper/blob/b58235d1cb42d5c7fd54122dc6b9f5e6c4110a75/src/TRACERxHelperFunctions.R#L1030  # nolint
                # let's just make sure we haven't created a subclonal mutation
                if (mut_cn / best_cn < 1L) {
                    if (best_cn > 1L) {
                        p <- prop_test_pvalues(
                            alt_count / 2L,
                            round(alt_count / (mut_cn / best_cn)),
                            0.5
                        )
                        if (p < 0.05) {
                            best_cn <- max(setdiff(best_cn, allCNs))
                        }
                    }
                }
            }
            c(out, best_cn)
        }, c(arg_list, list(alt_count = alt_counts, mut_cn = mutCopyNum)),
        MoreArgs = list(conipher = conipher)
    )
    data.table::transpose(out_list)
}

is_indel <- function(ref, alt) {
    nref <- data.table::fifelse(is.na(ref), 0L, nchar(ref))
    nalt <- data.table::fifelse(is.na(alt), 0L, nchar(alt))
    nref != nalt
}
