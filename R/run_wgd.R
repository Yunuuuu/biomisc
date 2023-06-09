#' Infer Whole-genome doubling
#'
#' @param seg_cnv A [data.frame][data.frame] obeject with following columns
#'  (arguments ends with "_field").
#' @param sample_field A string specifying the sample id column in seg_cnv.
#' @param contigs The Chromosome names to define Whole-genome doubling.
#' @param major_cn_field A string specifying the major_cn (allele-specific)
#'  column in seg_cnv. Default: "major_cn".
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
    major_cn_field = "major_cn", minor_cn_field = NULL,
    CNt_field = NULL, ploidy_field = "ploidy",
    chr_field = "chr", start_field = "startpos", end_field = "endpos",
    ref_cytoband = "hg38", arm_field = NULL, perm_times = 10000L) {
    assert_class(sample_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE,
        cross_msg = NULL
    )
    assert_class(major_cn_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = !is.null(minor_cn_field),
        cross_msg = NULL
    )
    assert_class(minor_cn_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = TRUE,
        cross_msg = NULL
    )
    assert_class(CNt_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = !(!is.null(minor_cn_field) && is.null(major_cn_field)),
        cross_msg = NULL
    )
    assert_class(ploidy_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        null_ok = is.null(minor_cn_field),
        cross_msg = NULL
    )
    assert_class(seg_cnv, is_class = "data.frame")
    if (is.null(minor_cn_field)) {
        group_fields <- sample_field
        other_fields <- major_cn_field
    } else {
        if (is.null(CNt_field)) {
            CNt_field <- "CNt"
            seg_cnv$CNt <- seg_cnv[[major_cn_field]] + seg_cnv[[minor_cn_field]]
        }
        group_fields <- c(sample_field, ploidy_field)
        other_fields <- c(minor_cn_field, CNt_field)
        assert_nest(seg_cnv, ploidy_field, group = sample_field)
    }
    out <- prepare_granges(
        data = seg_cnv,
        chr_field = chr_field,
        start_field = start_field,
        end_field = end_field,
        other_fields = c(group_fields, other_fields),
        keep.extra.columns = TRUE,
        ignore.strand = TRUE
    )
    assert_range_unique(out, group = sample_field)

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
    out <- map_seqnames(out, cytoband_seqstyle)
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
    out <- seg_to_arm(
        seg_cnv = out,
        arm_cytoband = arm_cytoband,
        arm_field = arm_field,
        group_fields = group_fields,
        other_fields = other_fields
    )
    if (is.null(minor_cn_field)) {
        data.table::setnames(out, major_cn_field, "major_cn")
        out[, width := as.double(width)]
        names(thresholds) <- paste0("Genome_prop_above_", thresholds)
        out <- out[, lapply(thresholds, function(threshold) {
            sum(width[major_cn >= threshold], na.rm = TRUE)
        }), by = c(sample_field, "chr", "arm", "arm_width")]
        out[, lapply(.SD, function(x) sum(x) / sum(arm_width)),
            by = sample_field, .SDcols = names(thresholds)
        ]
    } else {
        data.table::setnames(
            out, c(minor_cn_field, CNt_field, ploidy_field),
            c("arm_minor_ploidy", "arm_total_ploidy", "ploidy")
        )
        out <- out[!is.na(width),
            lapply(.SD, function(x) {
                matrixStats::weightedMedian(
                    x, w = width, na.rm = TRUE # styler: off
                )
            }),
            by = c(sample_field, "chr", "arm", "ploidy"),
            .SDcols = c("arm_minor_ploidy", "arm_total_ploidy")
        ]
        # nolint start
        out[, arm_major_ploidy := arm_total_ploidy - arm_minor_ploidy]
        out[, total_aber := sum(abs(arm_minor_ploidy - 1L), na.rm = TRUE) +
            sum(abs(arm_major_ploidy - 1L), na.rm = TRUE), by = sample_field]
        out[,
            c("A", "B") := lapply(.SD, function(x) {
                sign(x - 1L)
            }),
            .SDcols = c("arm_major_ploidy", "arm_minor_ploidy")
        ]
        out[total_aber > 0L,
            c("A_prob", "B_prob") := lapply(.SD, function(x) {
                abs(x - 1L) / total_aber
            }),
            .SDcols = c("A", "B")
        ]
        out <- out[,
            calculate_wgd(
                .SD, arm_major_ploidy, total_aber[[1L]], perm_times
            ),
            .SDcols = c("A", "B", "A_prob", "B_prob"),
            by = c(sample_field, "ploidy")
        ]
        out[, wGD := wgd_staus(pvalue, ploidy)][]
        # nolint end
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


# https://bitbucket.org/nmcgranahan/pancancerclonality/src/master/EstimateClonality/R/GD.functions.R
calculate_wgd <- function(events, arm_major_ploidy, total_aber, perm_times) {
    obs_wgd_event <- mean(arm_major_ploidy >= 2L)
    out_list <- list(total_aber = total_aber, wgd_event = obs_wgd_event)
    if (total_aber == 0L) {
        pvalue <- 1
    } else {
        if (obs_wgd_event == 1L) {
            pvalue <- 0
        } else {
            perm_wgd_events <- do.call(
                perm_wgd, c(
                    events,
                    list(total = total_aber, times = perm_times)
                )
            )
            pvalue <- mean(perm_wgd_events >= obs_wgd_event)
        }
        pvalue <- pvalue
    }
    c(pvalue = pvalue, out_list)
}

perm_wgd <- function(A, B, A_prob, B_prob, total, times = 1000L) {
    genotype <- cbind(A, B)
    vapply(seq_len(times), function(i, genotype) {
        idx <- sample(seq_along(c(A, B)), total,
            replace = TRUE, prob = c(A_prob, B_prob)
        )
        idx <- c(table(idx)) # name is index and value is the counts
        genotype[as.integer(names(idx))] <- 1L +
            genotype[as.integer(names(idx))] * idx
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

utils::globalVariables(c(
    "total_aber", "wGD", "pvalue", "ploidy", 
    "arm_major_ploidy", "arm_minor_ploidy",
    "arm_total_ploidy"
))
