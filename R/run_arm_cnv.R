#' Caculate Chromosome-arm-levels copy number variation
#'
#' @param seg_cnv A [data.frame][data.frame] obeject with CNV values in
#'   "cnv_field" column (absolute segment copy number determined by ABSOLUTE
#'   algorithm, which is the sum of allelic copy numbers, or relative segment
#'   copy number defined with -1 meaning del, 0 meaning neutral and 1 meaning
#'   amp).
#' @param sample_field A string indicates the sample Id column in seg_cnv.
#' @param cnv_field A scalar character gives the column containing CNV values.
#' @param chr_field,start_field,end_field A string specifying the column of the
#'  chromosome name, start positions and end positions of the genomic ranges in
#'  seg_cnv.
#' @param ref_cytoband A scalar string `"hg19"` or `"hg38"`, in this way, see
#'   [get_cytoband], or you can provided a self-defined
#'   [GenomicRanges][GenomicRanges::GenomicRanges] obeject containing the
#'   Cytoband reference.
#' @inheritParams get_arm_ranges
#' @param cnv_mode One of "rel" and "abs" correspongding to what Shukla, A and
#'   Cohen-Sharir have presented respectively.
#' @param ... Not used currently.
#' @param filter_centromere Whether to include or exclude segments across
#'   centromere, namely genomic ranges interseted with "acen" arm of
#'   ref_cytoband.  Default: `TRUE`.
#' @param threshold The fraction to define Chromosome arm-level aneuploidy
#'   profiling, Shukla, A. used 0.9 as the cut-off ("rel" cnv_mode).
#' @param ploidy The ploidy to define Chromosome arm-level aneuploidy profiling.
#'   Cohen-Sharir uses background ploidy ("abs" cnv_mode) derived from
#'   [run_absolute][ABSOLUTE::RunAbsolute] algorithm.
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
    seg_cnv, sample_field = NULL, cnv_field,
    chr_field = "chr", start_field = "startpos", end_field = "endpos",
    ref_cytoband = "hg38", arm_field = NULL, arms = c("p", "q"),
    cnv_mode = c("rel", "abs"),
    ..., filter_centromere = TRUE,
    ploidy = NULL, threshold = 0.9) {
    assert_pkg("GenomeInfoDb")
    assert_pkg("matrixStats")
    cnv_mode <- match.arg(cnv_mode)
    assert_class(sample_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL, null_ok = TRUE
    )
    assert_class(cnv_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    seg_cnv <- prepare_granges(
        data = seg_cnv,
        chr_field = chr_field,
        start_field = start_field,
        end_field = end_field,
        other_fields = c(sample_field, cnv_field),
        keep.extra.columns = TRUE,
        ignore.strand = TRUE
    )
    assert_range_unique(seg_cnv, group = sample_field)
    cnv_values <- S4Vectors::mcols(seg_cnv)[[cnv_field]]
    if (!is.numeric(cnv_values)) {
        cli::cli_abort("CNV values must be numeric")
    }

    if (cnv_mode == "abs") {
        if (!all(cnv_values >= 0L)) {
            cli::cli_abort("Copy number values must be positive or 0L in {.val {cnv_mode}} CNV mode")
        }
        ploidy <- ploidy %||% 2L
    } else if (cnv_mode == "rel") {
        S4Vectors::mcols(seg_cnv)[[cnv_field]] <- data.table::fcase(
            cnv_values < 0L, -1L, cnv_values > 0L, 1L,
            default = 0L
        )
    }
    if (rlang::is_scalar_character(ref_cytoband) &&
        any(ref_cytoband == c("hg19", "hg38"))) {
        ref_cytoband <- get_cytoband(ref_cytoband, add_arm = TRUE)
        arm_field <- "arm"
    } else if (!inherits(ref_cytoband, "GenomicRanges")) {
        cli::cli_abort(
            '{.arg ref_cytoband} must be a scalar character ("hg19" or "hg38"), or a self-defined {.cls GenomicRanges} object.'
        )
    }
    cytoband_seqstyle <- GenomeInfoDb::seqlevelsStyle(ref_cytoband)
    seg_cnv <- map_seqnames(seg_cnv, cytoband_seqstyle)

    arm_cytoband <- get_arm_ranges(ref_cytoband,
        arm_field = arm_field, arms = unique(c(arms, "acen"))
    )
    if (filter_centromere) {
        acen_region <- arm_cytoband[
            S4Vectors::mcols(arm_cytoband)[[arm_field]] == "acen"
        ]
        # Remove any segment that sligthly overlaps the centromere
        centromere_hits <- silent_expr(
            GenomicRanges::findOverlaps(seg_cnv, acen_region),
            "Each of the 2 combined objects has sequence levels not in the other:",
            fixed = TRUE
        )
        if (length(centromere_hits) > 0) {
            seg_cnv <- seg_cnv[-S4Vectors::queryHits(centromere_hits)]
        }
    }
    arm_cytoband <- arm_cytoband[
        S4Vectors::mcols(arm_cytoband)[[arm_field]] %in% arms
    ]
    out <- seg_to_arm(
        seg_cnv = seg_cnv,
        arm_cytoband = arm_cytoband,
        arm_field = arm_field,
        other_fields = c(sample_field, cnv_field)
    )
    data.table::setnames(out, cnv_field, "CNV")
    switch(cnv_mode,
        rel = out[, list(arm_cnv = as.integer(
            sum(abs(CNV) * width / arm_width, na.rm = TRUE) > threshold # nolint
        )), by = c(sample_field, "chr", "arm")],
        # https://github.com/broadinstitute/Aneuploidy_dependencies/blob/master/make_CCLE_arm_calls.R
        abs = out[,
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
    )
}

utils::globalVariables(c("CNV", "width", "arm_width"))
