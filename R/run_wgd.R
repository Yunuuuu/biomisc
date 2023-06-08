#' Infer Whole-genome doubling
#'
#' @param seg_cnv A [data.frame][data.frame] obeject with CNV values in
#'  "major_cn_field" .
#' @param sample_field A string specifying the sample id column in seg_cnv.
#' @param contigs The Chromosome names to define Whole-genome doubling.
#' @param major_cn_field A string specifying the major_cn (allele-specific)
#'  column in seg_cnv.
#' @param thresholds An integer vector specifying the thresholds to define
#'  Whole-genome doubling.
#' @param minor_cn_field,CNt_field,ploidy_field A string specifying the minor_cn
#'  (minor allele), total_cn (segmented absoltue copy number) and ploidy column
#'  in seg_cnv. Only used when `major_cn_field` is NULL.
#' @inheritParams run_arm_cnv
#' @param perm_times An integer specifying the times of simulation genome.
#' @references
#' - Bielski, C.M., Zehir, A., Penson, A.V. et al. Genome doubling shapes the
#'   evolution and prognosis of advanced cancers. Nat Genet 50, 1189–1195
#'   (2018).  <https://doi.org/10.1038/s41588-018-0165-1> 
#' - Tolerance of Whole-Genome Doubling Propagates Chromosomal Instability and
#'   Accelerates Cancer Genome Evolution.
#'   <https://doi.org/10.1158/2159-8290.CD-13-0285>
#' @export
run_wgd <- function(
    seg_cnv, sample_field = NULL,
    thresholds = 2:3, contigs = 1:22,
    major_cn_field = "major_cn",
    minor_cn_field = NULL, CNt_field = NULL, ploidy_field = NULL,
    chr_field = "chr", start_field = "startpos", end_field = "endpos",
    ref_cytoband = "hg38", arm_field = NULL, perm_times = 10000L) {
    assert_class(sample_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE,
        cross_msg = NULL
    )
    assert_class(major_cn_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE,
        cross_msg = NULL
    )
    assert_class(minor_cn_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = !is.null(major_cn_field),
        cross_msg = NULL
    )
    assert_class(CNt_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = !is.null(major_cn_field),
        cross_msg = NULL
    )
    assert_class(ploidy_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = !is.null(major_cn_field),
        cross_msg = NULL
    )

    seg_cnv <- prepare_granges(
        data = seg_cnv,
        chr_field = chr_field,
        start_field = start_field,
        end_field = end_field,
        other_fields = c(
            sample_field, major_cn_field,
            minor_cn_field, CNt_field, ploidy_field
        ),
        keep.extra.columns = TRUE,
        ignore.strand = TRUE
    )
    if (is.null(major_cn_field)) {
        assert_nest(seg_cnv, ploidy_field, group = sample_field)
    }
    assert_range_unique(seg_cnv, group = sample_field)

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
    contigs <- map_seqnames(as.character(contigs),
        cytoband_seqstyle,
        arg = "contigs"
    )
    arm_cytoband <- get_arm_ranges(ref_cytoband,
        arm_field = arm_field, arms = c("p", "q")
    )
    missing_contigs <- setdiff(contigs, GenomeInfoDb::seqnames(arm_cytoband)) # nolint
    if (length(missing_contigs) > 0L) {
        cli::cli_abort(c(
            "Cannot find all contigs in {.arg ref_cytoband}",
            x = "missing contig{?s}: {.val {missing_contigs}}"
        ))
    }
    arm_cytoband <- arm_cytoband[
        as.character(GenomeInfoDb::seqnames(arm_cytoband)) %chin% contigs
    ]
    seg_cnv <- seg_to_arm(
        seg_cnv = seg_cnv,
        arm_cytoband = arm_cytoband,
        arm_field = arm_field,
        other_fields = c(sample_field, major_cn_field, minor_cn_field)
    )
    if (!is.null(major_cn_field)) {
        data.table::setnames(seg_cnv, major_cn_field, "major_cn")
        seg_cnv <- seg_cnv[, c(
            arm_width = unique(arm_width),
            lapply(rlang::set_names(thresholds), function(threshold) {
                sum(width[major_cn >= threshold])
            })
        ), by = c(sample_field, "chr", "arm")]
        seg_cnv[, lapply(.SD, function(x) sum(x) / sum(arm_width)),
            by = sample_field, .SDcols = setdiff(
                names(seg_cnv), c(sample_field, "chr", "arm", "arm_width")
            )
        ]
        # seg_cnv[,
        #     perm_wgd(
        #         len = length(contigs), .SD, thresholds = thresholds,
        #         genome_width = genome_width, times = perm_times
        #     ),
        #     .SDcols = setdiff(names(seg_cnv), c(sample_field, "chr"))
        # ]
    } else {
        data.table::setnames(
            seg_cnv, c(minor_cn_field, CNt_field, ploidy_field),
            c("arm_minor_ploidy", "arm_total_ploidy", "ploidy")
        )
        out <- seg_cnv[,
            lapply(.SD, function(x) {
                matrixStats::weightedMedian(x, w = width, na.rm = TRUE)
            }),
            by = c(sample_field, "chr", "arm"),
            .SDcols = c("arm_minor_ploidy", "arm_total_ploidy")
        ]
        out[, arm_major_ploidy := arm_total_ploidy - arm_minor_ploidy] # nolint
        out[, total_aber := sum(abs(arm_minor_ploidy - 1L)) + # nolint
            sum(abs(arm_major_ploidy - 1L)), by = sample_field] # nolint
        out[,
            c("A", "B") := lapply(.SD, function(x) {
                sign(x - 1L)
            }),
            .SDcols = c("arm_major_ploidy", "arm_minor_ploidy")
        ]
        out[,
            c("A_prob", "B_prob") := lapply(.SD, function(x) {
                abs(x) / total_aber # nolint
            }),
            .SDcols = c("A", "B")
        ]
        out <- out[,
            {
                obs_wgd_event <- mean(arm_major_ploidy >= 2L) # nolint
                if (obs_wgd_event == 1L) {
                    pvalue <- 0L
                } else {
                    perm_wgd_events <- do.call(
                        perm_wgd, c(.SD, list(
                            total = total_aber[[1L]], times = perm_times # nolint
                        ))
                    )
                    pvalue <- mean(perm_wgd_events >= obs_wgd_event)
                }
                list(pvalue = pvalue, ploidy = ploidy[[1L]]) # nolint
            },
            .SDcols = c("A", "B", "A_prob", "B_prob"),
            by = sample_field
        ]
        out[, wGD := wgd_staus(pvalue, ploidy)][] # nolint
    }
}

# we shoudn't permutation in a single sample but should in a cohort.
# perm_wgd <- function(size, events, thresholds, genome_width, times = 1000L) {
#     totals <- nrow(events)
#     perm_stats_list <- lapply(seq_len(times), function(i) {
#         idx <- sample.int(n = totals, size = size, replace = TRUE)
#         colSums(events[idx]) / genome_width
#     })
#     perm_stats_list <- data.table::transpose(perm_stats_list)
#     pvalues <- mapply(function(stat, perm_stats) {
#         mean(perm_stats >= stat)
#     }, stat = actual_stats, perm_stats = perm_stats_list, USE.NAMES = FALSE)
#     list(thretholds = thresholds, wgd_prob = actual_stats, pvalues = pvalues)
# }

perm_wgd <- function(A, B, A_prob, B_prob, total, times = 1000L) {
    genotype <- cbind(A, B)
    vapply(seq_len(times), function(i, genotype) {
        idx <- sample(seq_along(c(A, B)), total,
            replace = TRUE, prob = c(A_prob, B_prob)
        )
        idx <- c(table(idx)) # name is index and value is the counts
        genotype[as.integer(names(idx))] <- 1L + genotype * idx
        genome <- pmax(genotype[, 1L, drop = TRUE], genotype[, 2L, drop = TRUE])
        mean(genome >= 2L)
    }, numeric(1L), genotype = genotype)
}

wgd_staus <- function(pvalue, ploidy) {
    ploidy <- round(ploidy)
    data.table::fcase(
        ploidy >= 6L & pvalue <= 1L, "GD",
        ploidy == 5L & pvalue <= 0.5, "GD",
        ploidy == 4L & pvalue <= 0.05, "GD",
        ploidy <= 3L & pvalue <= 0.001, "GD",
        default = "nGD"
    )
}
