#' Estimating CCF (absolute and phylogenetic)
#'
#' @param mut_cn_data A data.frame with mutation and copy number data. Copy
#'  number often contain subclonal copy number as described in
#'  [CONIPHER](https://github.com/McGranahanLab/CONIPHER-wrapper/blob/b58235d1cb42d5c7fd54122dc6b9f5e6c4110a75/src/TRACERxHelperFunctions.R#L1).
#' @param patient_field A string specifying the patient column. If NULL, all
#'  data will be regarded from the same patient.
#' @param sample_field A string specifying the sample column. If NULL, all data
#'  will be regarded from the same sample. This is used to confirm every sample
#'  have the same `purity` or `gender`.
#' @seealso [run_ccf]
#' @note Just for internal usage.
#' @noRd
estimate_phylo_ccf <- function(
    mut_cn_data, patient_field = NULL, sample_field = NULL,
    chr_field = NULL, pos_field = NULL, ref_field = NULL,
    alt_field = NULL, indel_field = NULL,
    min_vaf_to_explain = NULL, conipher = FALSE) {
    assert_(
        ref_field, rlang::is_scalar_character,
        "a string",
        null_ok = is.null(indel_field), show_length = TRUE
    )
    assert_(
        alt_field, rlang::is_scalar_character,
        "a string",
        null_ok = TRUE, show_length = TRUE
    )
    assert_(
        indel_field, rlang::is_scalar_character,
        "a string",
        null_ok = TRUE, show_length = TRUE
    )
    assert_(min_vaf_to_explain, is_scalar_numeric,
        "a number",
        null_ok = TRUE, show_length = TRUE
    )
    assert_(conipher, rlang::is_scalar_logical,
        "a scalar {.cls logical}",
        null_ok = TRUE, show_length = TRUE
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

    # nolint start
    out[, expVAF := calculate_vaf(
        1L, purity, major_raw + minor_raw, normal_cn
    )]
    out[, obsVAF := alt_counts / (alt_counts + ref_counts)]

    # following estimate phyloCCF and absoluteCCF
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
    if (conipher) {
        out[
            ,
            ..tmp_matched_rows.. := ..tmp_matched_rows.. |
                ..tmp_best_cn.. >= (1 - 0.05)
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
    out[
        ..tmp_amp_mut_pvalue.. <= 0.05 & mutCopyNum > 1L,
        c("phyloCCF", "phyloCCF_lower", "phyloCCF_higher", "expVAF") := list(
            mutCopyNum / ..tmp_best_cn..,
            phyloCCF_lower / ..tmp_best_cn..,
            phyloCCF_higher / ..tmp_best_cn..,
            expVAF * ..tmp_best_cn..
        )
    ]
    out[, ..tmp_amp_mut_pvalue.. := NULL]
    out[, ..tmp_best_cn.. := NULL]
    # nolint end

    # finally, let's sort out 'missing' ones
    out[
        alt_counts == 0L, # nolint
        c("no.chrs.bearing.mut", "expVAF", "absCCF") := list(0L, 0L, 0L)
    ]

    # correct Indel CCF
    # to ensure potentially unreliable VAFs of indels did not lead to
    # separate mutation clusters, each estimated indel CCF was multiplied by a
    # region specific correction factor.  Assuming the majority of ubiquitous
    # mutations, present in all regions, are clonal, the region specific
    # correction factor was calculated by dividing the median mutation CCF of
    # ubiquitous mutations by the median indel CCF of ubiquitous indels.
    # SNV: width(ref) == width(alt) including snv, dbs, mbs
    if (is.null(indel_field)) {
        if (!is.null(ref_field) && !is.null(alt_field)) {
            out$..tmp_is_indel.. <- is_indel(out[[ref_field]], out[[alt_field]])
        } else if (!is.null(ref_field) || !is.null(alt_field)) {
            cli::cli_abort(
                "both {.arg ref_field} and {.arg alt_field} must be specified to correct indel phyloCCF, Otherwise, you should leave both as NULL"
            )
        }
    } else {
        if (!is.logical(out[[indel_field]])) {
            cli::cli_abort("{indel_field} column of {.arg mut_data} must be logical")
        }
        data.table::setnames(out, indel_field, "..tmp_is_indel..")
    }
    if (!is.null(out$..tmp_is_indel..) && any(out$..tmp_is_indel..)) {
        # we keep backup sample numbers of every patient
        if (is.null(sample_field)) {
            out[, ..tmp_nsamples.. := 1L] # nolint
        } else {
            if (is.null(patient_field)) {
                out[, ..tmp_nsamples.. := length(unique(.SD[[1L]])), # nolint
                    .SDcols = sample_field
                ]
            } else {
                out[, ..tmp_nsamples.. := length(unique(.SD[[1L]])), # nolint
                    .SDcols = sample_field, by = patient_field
                ]
            }
        }
        # then, we identify ubiquitous mutations
        ubiq_mut <- out[, .SD[sum(obsVAF > 5L) == ..tmp_nsamples..[1L]], # nolint
            by = c(patient_field, chr_field, pos_field, ref_field)
        ][!is.na(phyloCCF)] # nolint

        if (nrow(ubiq_mut)) {
            cli::cli_inform("Correcting Indel phyloCCF")
            # nolint start
            # for each sample, we calcualte indel correction factor based on
            # ubiquitous mutation
            if (is.null(sample_field) && is.null(patient_field)) {
                ..tmp_indel_cf.. <- ubiq_mut[
                    ,
                    stats::median(phyloCCF[!..tmp_is_indel..]) /
                        stats::median(phyloCCF[..tmp_is_indel..])
                ]
                out[, ..tmp_correction_factor.. := ..tmp_indel_cf..]
            } else {
                ..tmp_indel_cf.. <- ubiq_mut[, list(
                    ..tmp_correction_factor.. =
                        stats::median(phyloCCF[!..tmp_is_indel..]) /
                            stats::median(phyloCCF[..tmp_is_indel..])
                ), by = c(patient_field, sample_field)]
                out <- ..tmp_indel_cf..[out,
                    on = c(patient_field, sample_field)
                ]
            }
            # for each indel, we use correction factor to
            # correct phyloCCF and mutCopyNum
            out[
                (..tmp_is_indel..),
                phyloCCF := phyloCCF * ..tmp_correction_factor..
            ]
            out[
                (..tmp_is_indel..),
                mutCopyNum := mutCopyNum * ..tmp_correction_factor..
            ]
            out[, ..tmp_correction_factor.. := NULL]
        }
        out[, ..tmp_is_indel.. := NULL]
        out[, ..tmp_nsamples.. := NULL]
        # nolint end
    }
    out
}
utils::globalVariables(c(
    "sample_id", "expVAF", "obsVAF", "purity",
    "fracA", "fracB", "fracC", "fracD",
    "fracMaj1", "fracMaj2", "fracMin1", "fracMin2",
    "major_raw", "minor_raw", "mutCopyNum", "nAraw", "nBraw",
    "nMajor", "nMinor", "minor_cn",
    "nMaj1", "nMaj2", "nMaj_A", "nMaj_B",
    "nMaj_C", "nMaj_D", "nMin1", "nMin2",
    "nMin_A", "nMin_B", "nMin_C", "nMin_D",
    "normal_cn", "CNn", "CNt", "phyloCCF",
    "phyloCCF_higher", "phyloCCF_lower", "absCCF_higher",
    "ref_counts", "startpos", "alt_counts", "whichFrac",
    "..tmp_operated_rows..", "..tmp_matched_rows..",
    "..tmp_mut_multi..", "..tmp_is_indel..",
    "..tmp_is_subclone..", "..tmp_amp_mut_pvalue..",
    "..tmp_explained_by_cn_pvalue..", "..tmp_expProp..",
    "..tmp_best_cn..", "..tmp_correction_factor..", "..tmp_nsamples.."
))

#' min.subclonal was set to 0.1 in
#' https://bitbucket.org/nmcgranahan/clonalneoantigenanalysispipeline/src/master/
#' @return A data.table
#' @keywords validated
#' @noRd
define_subclone_cn <- function(seg, min_subclonal = 0.01) {
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
    if (length(counts) == 0L) {
        return(rep_len(list(numeric()), 3L))
    }
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

is_indel <- function(ref, alt) {
    nref <- data.table::fifelse(is.na(ref), 0L, nchar(ref))
    nalt <- data.table::fifelse(is.na(alt), 0L, nchar(alt))
    nref != nalt
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
