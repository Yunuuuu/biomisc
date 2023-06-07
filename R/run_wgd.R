#' Infer Whole-genome doubling
#'
#' @param seg_cnv A [data.frame][data.frame] obeject with CNV values in
#'  "major_cn_field" .
#' @param sample_field A string specifying the sample id column in seg_cnv.
#' @param thresholds An integer vector specifying the thresholds to define
#'  Whole-genome doubling.
#' @param contigs The Chromosome names to define Whole-genome doubling.
#' @param major_cn_field A string specifying the major_cn (allele-specific)
#'  column in seg_cnv.
#' @inheritParams run_arm_cnv
#' @references Bielski, C.M., Zehir, A., Penson, A.V. et al. Genome doubling
#' shapes the evolution and prognosis of advanced cancers. Nat Genet 50,
#' 1189–1195 (2018). https://doi.org/10.1038/s41588-018-0165-1
#' @export
run_wgd <- function(
    seg_cnv, sample_field = NULL,
    thresholds = 2:3, contigs = 1:22,
    major_cn_field = "major_cn",
    chr_field = "chr", start_field = "startpos", end_field = "endpos",
    ref_cytoband = "hg38", arm_field = NULL) {
    assert_pkg("GenomicRanges")
    assert_pkg("GenomeInfoDb")
    assert_class(sample_field, rlang::is_scalar_character,
        "scalar {.cls character}", null_ok = TRUE,
        cross_msg = NULL
    )
    assert_class(major_cn_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    assert_class(chr_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    assert_class(start_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    assert_class(end_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        cross_msg = NULL
    )
    assert_df_with_columns(seg_cnv, c(
        sample_field, major_cn_field, chr_field, start_field, end_field
    ))
    # prepare seg CNV data --------------------------
    seg_cnv <- GenomicRanges::makeGRangesFromDataFrame(
        seg_cnv,
        keep.extra.columns = TRUE,
        seqnames.field = chr_field,
        start.field = start_field,
        end.field = end_field,
        ignore.strand = TRUE
    )
    # prepare cytoband data ------------------------
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
    contigs <- map_seqnames(as.character(contigs), cytoband_seqstyle)
    arm_cytoband <- get_arm_ranges(ref_cytoband,
        arm_field = arm_field, arms = c("p", "q")
    )
    missing_contigs <- setdiff(contigs, GenomeInfoDb::seqnames(arm_cytoband)) # nolint
    if (length(missing_contigs) > 0L) {
        cli::cli_abort(c(
            "Cannot find all contigs in {.arg ref_cytoband}",
            x = "missing contigs: {.val {missing_contigs}}"
        ))
    }
    arm_cytoband <- arm_cytoband[
        as.character(GenomeInfoDb::seqnames(arm_cytoband)) %chin% contigs
    ]
    genome_width <- sum(GenomicRanges::width(arm_cytoband))
    out <- seg_to_arm(
        seg_cnv = seg_cnv,
        arm_cytoband = arm_cytoband,
        arm_field = arm_field,
        other_fields = c(sample_field, major_cn_field)
    )
    if (is.null(sample_field)) {
        out <- out[contigs, on = "chr"]
    } else {
        out <- out[, .SD[contigs, on = "chr"], by = c(sample_field)]
    }
    data.table::setnames(out, major_cn_field, "major_cn")
    out <- out[, lapply(thresholds, function(threshold) {
        sum(width[major_cn >= threshold])
    }), by = c(sample_field, "chr")]
    out[,
        perm_wgd(
            chr = chr, .SD, thresholds = thresholds,
            genome_width = genome_width
        ),
        .SDcols = setdiff(names(out), c(sample_field, "chr")),
        by = sample_field
    ]
}

perm_wgd <- function(chr, events, thresholds, genome_width, times = 1000L) {
    len <- length(chr)
    actual_stats <- colSums(events) / genome_width
    perm_stats_list <- lapply(seq_len(times), function(i) {
        idx <- sample.int(len, replace = TRUE)
        colSums(events[idx]) / genome_width
    })
    perm_stats_list <- data.table::transpose(perm_stats_list)
    pvalues <- mapply(function(stat, perm_stats) {
        1L - stats::ecdf(perm_stats)(stat)
    }, stat = actual_stats, perm_stats = perm_stats_list, USE.NAMES = FALSE)
    list(thretholds = thresholds, wgd_events = actual_stats, pvalues = pvalues)
}
