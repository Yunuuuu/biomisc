#' Caculate Chromosome-arm-levels copy number variation
#'
#' @param seg_data A [data.frame][data.frame] obeject with segmented chromosome
#' data.
#' @param sample_field A string indicates the sample Id column in seg_data.
#' @param other_fields A character vector indicates other column in seg_data to
#' keep in the final result.
#' @param chr_field,start_field,end_field A string specifying the column of the
#'  chromosome name, start positions and end positions of the genomic ranges in
#'  seg_data.
#' @param ref_cytoband A scalar string `"hg19"` or `"hg38"`, in this way, see
#'   [get_cytoband], or you can provided a self-defined
#'   [GenomicRanges][GenomicRanges::GenomicRanges] obeject containing the
#'   Cytoband reference.
#' @param contigs A character vector specifying the chromosome to summarize.
#' @inheritParams get_arm_ranges
#' @param ... Not used currently.
#' @param filter_centromere Whether to include or exclude segments across
#'   centromere, namely genomic ranges interseted with "acen" arm of
#'   ref_cytoband.  Default: `TRUE`.
#' @param pruning_mode When some of the seqlevels to drop from `ref_cytoband`
#'  are in use (i.e. have ranges on them), the ranges on these sequences need to
#'  be removed before the seqlevels can be dropped. We call this pruning. The
#'  pruning_mode argument controls how to prune x. See
#'  [seqinfo][GenomeInfoDb::seqinfo] pruning.mode
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @return A [data.table][data.table::data.table] containing Chromosome-arm
#'  level data.
#' @export
summarize_arm <- function(
    seg_data, sample_field = NULL, other_fields = NULL,
    chr_field = "chr", start_field = "startpos", end_field = "endpos",
    ref_cytoband = "hg38", contigs = NULL,
    arm_field = NULL, arms = c("p", "q"), ...,
    filter_centromere = TRUE, pruning_mode = "error") {
    assert_pkg("GenomeInfoDb")
    assert_pkg("S4Vectors")
    assert_class(sample_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL, null_ok = TRUE
    )
    assert_class(filter_centromere, rlang::is_scalar_logical,
        "scalar {.cls logical}",
        cross_msg = NULL, null_ok = TRUE
    )
    seg_ranges <- prepare_granges(
        data = seg_data,
        chr_field = chr_field,
        start_field = start_field,
        end_field = end_field,
        other_fields = c(other_fields, sample_field),
        keep.extra.columns = TRUE,
        ignore.strand = TRUE
    )
    assert_range_unique(seg_ranges, group = sample_field)
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
    if (is.null(contigs)) {
        contigs <- GenomeInfoDb::seqnames(ref_cytoband)
    } else {
        contigs <- as.character(contigs)
        contigs <- map_seqnames(contigs, cytoband_seqstyle,
            style_arg = "ref_cytoband"
        )
        ref_cytoband <- GenomeInfoDb::keepSeqlevels(
            ref_cytoband, contigs,
            pruning.mode = pruning_mode
        )
    }
    seg_ranges <- map_seqnames(seg_ranges, cytoband_seqstyle,
        arg = "seg_data", style_arg = "ref_cytoband"
    )

    arm_cytoband <- get_arm_ranges(ref_cytoband,
        arm_field = arm_field, arms = unique(c(arms, "acen"))
    )
    if (filter_centromere) {
        acen_region <- arm_cytoband[
            S4Vectors::mcols(arm_cytoband)[[arm_field]] == "acen"
        ]
        # Remove any segment that sligthly overlaps the centromere
        centromere_hits <- silent_expr(
            GenomicRanges::findOverlaps(seg_ranges, acen_region),
            "Each of the 2 combined objects has sequence levels not in the other:",
            fixed = TRUE
        )
        if (length(centromere_hits) > 0) {
            seg_ranges <- seg_ranges[-S4Vectors::queryHits(centromere_hits)]
        }
    }
    arm_cytoband <- arm_cytoband[
        S4Vectors::mcols(arm_cytoband)[[arm_field]] %in% arms
    ]
    out <- seg_to_arm_seg(
        seg_ranges = seg_ranges,
        arm_cytoband = arm_cytoband,
        arm_field = arm_field,
        group_fields = sample_field,
        other_fields = other_fields,
        arg1 = "seg_data", arg2 = "ref_cytoband"
    )
    out[]
}

seg_to_arm_seg <- function(seg_ranges, arm_cytoband, arm_field, group_fields = NULL, other_fields = character(), arg1 = rlang::caller_arg(seg_ranges), arg2 = rlang::caller_arg(arm_cytoband), call = parent.frame()) {
    # find ouverlap index --------------------------------
    overlap_hits <- silent_expr(
        GenomicRanges::findOverlaps(
            seg_ranges, arm_cytoband,
            type = "any"
        ),
        "Each of the 2 combined objects has sequence levels not in the other:",
        fixed = TRUE
    )
    if (length(overlap_hits) == 0L) {
        cli::cli_abort("No match between {.arg {arg1}} and {.arg {arg2}}",
            call = call
        )
    }
    seg_hits <- S4Vectors::queryHits(overlap_hits)
    cytoband_hits <- S4Vectors::subjectHits(overlap_hits)

    # extract intersected ranges ------------------------
    intersect_region <- silent_expr(
        GenomicRanges::pintersect(
            seg_ranges[seg_hits], arm_cytoband[cytoband_hits],
            drop.nohit.ranges = FALSE,
            ignore.strand = FALSE,
            strict.strand = FALSE
        ),
        "Each of the 2 combined objects has sequence levels not in the other:",
        fixed = TRUE
    )

    # define arm-level cnv -------------------------------
    arm_ranges <- arm_cytoband[cytoband_hits]
    out <- data.table::data.table(
        chr = as.factor(GenomicRanges::seqnames(intersect_region)),
        arm = S4Vectors::mcols(arm_ranges)[[arm_field]],
        width = GenomicRanges::width(intersect_region)
    )

    # ensure everything in the output
    other_fields <- c(other_fields, group_fields)
    if (length(other_fields)) {
        other_data <- data.table::as.data.table(seg_ranges[seg_hits])
        other_data <- other_data[, .SD, .SDcols = other_fields]
    }
    out <- cbind(out, other_data)

    # ensure every reference arm is in the out result
    ref_cytoband_dt <- unique(data.table::as.data.table(arm_cytoband))
    data.table::setnames(
        ref_cytoband_dt,
        c("seqnames", arm_field, "width"), c("chr", "arm", "arm_width")
    )
    ref_cytoband_dt <- ref_cytoband_dt[, c("chr", "arm", "arm_width")]
    if (is.null(group_fields)) {
        out <- out[ref_cytoband_dt,
            on = c("chr", "arm"),
            allow.cartesian = FALSE
        ]
    } else {
        out <- out[,
            .SD[ref_cytoband_dt, on = c("chr", "arm"), allow.cartesian = FALSE],
            by = group_fields
        ]
    }
    # finally, move arm_width next to arm column
    data.table::setcolorder(out, "arm_width", after = "arm")
}
