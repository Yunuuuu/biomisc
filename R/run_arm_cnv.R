#' Caculate Chromosome-arm-levels copy number variation
#'
#' @inheritParams summarize_arm
#' @param seg_data A [data.frame][data.frame] obeject with CNV values in
#'   "cnv_field" column (absolute segment copy number determined by ABSOLUTE
#'   algorithm, which is the sum of allelic copy numbers, or relative segment
#'   copy number defined with -1 meaning del, 0 meaning neutral and 1 meaning
#'   amp).
#' @param cnv_field A scalar character gives the column containing CNV values.
#' @param cnv_mode One of "rel" and "abs" correspongding to what Shukla, A and
#'   Cohen-Sharir have presented respectively.
#' @param threshold The fraction to define Chromosome arm-level aneuploidy
#'   profiling, Shukla, A. used 0.9 as the cut-off ("rel" cnv_mode).
#' @param ploidy_field A single numeric indicates the ploidy to define
#'   Chromosome arm-level aneuploidy profiling.  Cohen-Sharir uses background
#'   ploidy ("abs" cnv_mode) derived from [run_absolute][ABSOLUTE::RunAbsolute]
#'   algorithm. Or you can also provide a string to specifying the ploidy column
#'   in `seg_cnv`.
#' @inheritDotParams summarize_arm -seg_data -sample_field -other_fields -group_fields
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @return A [data.table][data.table::data.table] containing
#' Chromosome-arm-levels copy number.
#' @references  \itemize{
#'  \item Cohen-Sharir, Y., McFarland, J.M., Abdusamad, M. et al. Aneuploidy
#'    renders cancer cells vulnerable to mitotic checkpoint inhibition. Nature
#'    590, 486–491 (2021).  \url{https://doi.org/10.1038/s41586-020-03114-6}
#'  \item Shukla, A., Nguyen, T.H.M., Moka, S.B. et al. Chromosome arm
#'    aneuploidies shape tumour evolution and drug response. Nat Commun 11, 449
#'    (2020).  \url{https://doi.org/10.1038/s41467-020-14286-0}
#' }
#' @export
run_arm_cnv <- function(
    seg_data, sample_field = NULL, cnv_field, cnv_mode = c("rel", "abs"),
    ploidy_field = 2L, threshold = 0.9, ...) {
    assert_pkg("GenomeInfoDb")
    assert_pkg("matrixStats")
    cnv_mode <- match.arg(cnv_mode)
    assert_(sample_field, rlang::is_scalar_character,
        "a string",
        null_ok = TRUE, show_length = TRUE
    )
    assert_(
        cnv_field, rlang::is_scalar_character,
        "a string",
    )
    group_fields <- NULL
    if (cnv_mode == "abs") {
        assert_(ploidy_field,
            function(x) {
                rlang::is_scalar_character(x) || is_scalar_numeric(x)
            }, "a number or a string",
            null_ok = TRUE, show_length = TRUE
        )
        if (is.character(ploidy_field)) {
            assert_data_frame_hierarchy(seg_data, ploidy_field, sample_field)
            group_fields <- ploidy_field
        }
    } else {
        assert_(
            threshold,
            is_scalar_numeric, "a number",
            show_length = TRUE
        )
    }
    out <- summarize_arm(
        seg_data = seg_data,
        sample_field = sample_field,
        other_fields = cnv_field,
        group_fields = group_fields,
        ...
    )
    data.table::setnames(out, cnv_field, "CNV")
    if (!is.numeric(out$CNV)) {
        cli::cli_abort("CNV values must be numeric")
    }
    if (cnv_mode == "abs") {
        if (!all(out$CNV >= 0L, na.rm = TRUE)) {
            cli::cli_abort("Copy number values must be positive or 0L in {.val {cnv_mode}} CNV mode")
        }
    } else {
        out$CNV <- data.table::fcase(
            out$CNV < 0L, -1L, out$CNV > 0L, 1L,
            default = 0L
        )
    }
    if (cnv_mode == "abs") {
        if (is.character(ploidy_field)) {
            data.table::setnames(out, ploidy_field, "ploidy")
        } else {
            out[, ploidy := ploidy_field]
        }
        # https://github.com/broadinstitute/Aneuploidy_dependencies/blob/master/make_CCLE_arm_calls.R
        out[,
            {
                median_weighted <- matrixStats::weightedMedian(
                    CNV, w = width, na.rm = TRUE # nolint styler: off
                )
                median_weighted <- round(median_weighted, digits = 0L)
                ploidy <- round(ploidy, digits = 0L)
                list(arm_cnv = as.integer(sign(median_weighted - ploidy)))
            },
            by = c(sample_field, "chr", "arm")
        ]
    } else {
        out[,
            list(
                arm_cnv = as.integer(
                    sum(abs(CNV) * width / arm_width, na.rm = TRUE) > threshold
                )
            ),
            by = c(sample_field, "chr", "arm")
        ]
    }
}
