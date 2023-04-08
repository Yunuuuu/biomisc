#' Caculate Chromosome-arm-levels copy number variation
#'
#' @param seg_cnv A [GenomicRanges][GenomicRanges::GenomicRanges] obeject with
#'   CNV values in "cnv_col" column of metadata (absolute segment copy number
#'   determined by ABSOLUTE algorithm, which is the sum of allelic copy numbers,
#'   or relative segment copy number defined with -1 meaning del, 0 meaning
#'   neutral and 1 meaning amp) and samples IDs in "sample_id_col" column.
#' @param cnv_col A scalar character gives the column containing CNV values.
#' @param ref_cytoband A [GenomicRanges][GenomicRanges::GenomicRanges] obeject
#'   containing the Cytoband reference, It can be a scalar character `"hg19"` or
#'   `"hg38"`, in this way, see [get_cytoband], or you can provided a
#'   self-defined [GenomicRanges][GenomicRanges::GenomicRanges] obeject.
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
#' @return A data.frame containing Chromosome-arm-levels copy number.
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
    seg_cnv, cnv_col,
    ref_cytoband = "hg38", arm_col = NULL,
    cnv_mode = c("rel", "abs"),
    ..., filter_centromere = TRUE,
    ploidy = NULL,
    threshold = 0.9) {
    assert_pkg("S4Vectors")
    assert_pkg("GenomicRanges")
    assert_pkg("GenomeInfoDb")
    assert_pkg("matrixStats")
    assert_class(seg_cnv, "GenomicRanges")

    cnv_mode <- match.arg(cnv_mode)

    if (!rlang::is_scalar_character(cnv_col)) {
        cli::cli_abort("{.arg cnv_col} must be a scalar character")
    } else if (!any(cnv_col == colnames(S4Vectors::mcols(seg_cnv)))) {
        cli::cli_abort("Cannot find {.val {cnv_col}} in the metadata column of {.arg seg_cnv}")
    }
    cnv_values <- S4Vectors::mcols(seg_cnv)[[cnv_col]]
    if (!is.numeric(cnv_values)) {
        cli::cli_abort("CNV values must be numeric")
    }

    if (cnv_mode == "abs") {
        if (!all(cnv_values >= 0L)) {
            cli::cli_abort("Copy number values must be positive or 0L in {.val {cnv_mode}} CNV mode")
        }
        if (is.null(ploidy)) {
            ploidy <- 2L
        }
    } else if (cnv_mode == "rel") {
        S4Vectors::mcols(seg_cnv)[[cnv_col]] <- data.table::fcase(
            cnv_values < 0L, -1L, cnv_values > 0L, 1L,
            default = 0L
        )
    }

    if (all(GenomeInfoDb::seqlevelsStyle(seg_cnv) != "UCSC")) {
        cli::cli_inform("Mapping seqnames of {.arg seg_cnv} to UCSC style")
        GenomeInfoDb::seqlevels(seg_cnv) <- "UCSC"
    }

    if (rlang::is_scalar_character(ref_cytoband) &&
        any(ref_cytoband == c("hg19", "hg38"))) {
        ref_cytoband <- get_cytoband(ref_cytoband, add_arm = TRUE)
        arm_col <- "arm"
    } else if (!inherits(ref_cytoband, "GenomicRanges")) {
        cli::cli_abort(
            '{.arg ref_cytoband} must be a scalar character ("hg19" or "hg38"), or a self-defined {.cls GenomicRanges} object.'
        )
    }
    arm_cytoband <- get_arm_ranges(ref_cytoband, arm_col = arm_col)
    if (filter_centromere) {
        acen_region <- arm_cytoband[arm_cytoband$arm == "acen"]
        # Remove any segment that sligthly overlaps the centromere
        centromere_hits <- GenomicRanges::findOverlaps(seg_cnv, acen_region)
        if (length(centromere_hits) > 0) {
            seg_cnv <- seg_cnv[-S4Vectors::queryHits(centromere_hits)]
        }
    }
    arm_cytoband <- arm_cytoband[arm_cytoband$arm %chin% c("p", "q")]
    seg_to_arm_cnv(
        seg_cnv = seg_cnv,
        cnv_col = cnv_col,
        arm_cytoband = arm_cytoband,
        cnv_mode = cnv_mode,
        threshold = threshold,
        ploidy = ploidy
    )
}

seg_to_arm_cnv <- function(seg_cnv, cnv_col, arm_cytoband, cnv_mode = "rel", ..., threshold, ploidy) {
    # find ouverlap index --------------------------------
    overlap_hits <- GenomicRanges::findOverlaps(
        seg_cnv, arm_cytoband,
        type = "any"
    )
    seg_hits <- S4Vectors::queryHits(overlap_hits)
    cytoband_hits <- S4Vectors::subjectHits(overlap_hits)

    # extract intersected ranges ------------------------
    intersect_region <- GenomicRanges::pintersect(
        seg_cnv[seg_hits], arm_cytoband[cytoband_hits],
        drop.nohit.ranges = FALSE,
        ignore.strand = FALSE,
        strict.strand = FALSE
    )

    # define arm-level cnv -------------------------------
    arm_ranges <- arm_cytoband[cytoband_hits]
    out <- data.table::data.table(
        seqnames = as.factor(GenomicRanges::seqnames(intersect_region)),
        width = GenomicRanges::width(intersect_region),
        # though we need the abolsute copy number only in rel cnv_mode
        # since the copy numbers of abs cnv_mode are always larger than zero
        # it won't hurt to do it for both cnv_mode.
        CNV = abs(S4Vectors::mcols(intersect_region)[[cnv_col]]),
        arm = S4Vectors::mcols(arm_ranges)[["arm"]],
        arm_width = GenomicRanges::width(arm_ranges)
    )
    out <- switch(cnv_mode,
        rel = out[, list(arm_cnv = as.integer(
            sum(CNV * width / arm_width, na.rm = TRUE) > threshold # nolint
        )), by = c("seqnames", "arm")],
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
            by = c("seqnames", "arm")
        ]
    )
    data.table::setDF(out)
}

utils::globalVariables(c(
    "CNV", "width", "arm_width"
))
