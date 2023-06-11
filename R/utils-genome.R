#' Get arm-level ranges
#' @param ref_ranges A [`GenomicRanges`][GenomicRanges::GRanges-class] object to
#'   combine into arm-level ranges. It's easy to use [get_cytoband] to get such
#'   Genomic reference ranges.
#' @param arm_field A scalar string indicates the column containing the
#'   chromosome arm. If `NULL`, the internal will search the first column starts
#'   with "arm" (ignore letter case). Only values of "p", "q", "acen" and "" (no
#'   arm) are supported.
#' @param arms A character specifying the arms to return.
#' @return A [GenomicRanges][GenomicRanges::GRanges-class] object containing
#' arm-level ranges.
#' @seealso [get_cytoband]
#' @export
get_arm_ranges <- function(ref_ranges, arm_field = NULL, arms = c("p", "q")) {
    assert_pkg("S4Vectors")
    assert_pkg("GenomicRanges")
    assert_pkg("GenomeInfoDb")
    assert_class(ref_ranges, "GenomicRanges")

    if (is.null(arm_field)) {
        arm_field <- grep("^arm", names(S4Vectors::mcols(ref_ranges)),
            ignore.case = TRUE, perl = TRUE,
            value = TRUE
        )
        if (!length(arm_field)) {
            cli::cli_abort(c(
                "Cannot determine right {.arg arm_field}",
                "i" = "try to set {.arg arm_field}"
            ))
        }
        arm_field <- arm_field[[1L]]
    } else if (!rlang::is_scalar_character(arm_field)) {
        cli::cli_abort("{.arg arm_field} should be a scalar string")
    }
    arm_values <- S4Vectors::mcols(ref_ranges)[[arm_field]]
    arm_levels <- c("p", "acen", "q", "")

    if (!all(as.character(arm_values) %chin% arm_levels)) {
        cli::cli_abort("Only values of {.val {arm_levels}} are supported in column specified in {.arg arm_field}.")
    }
    arms <- as.character(arms)
    if (!all(arms %chin% arm_levels)) {
        cli::cli_abort("Only values of {.val {arm_levels}} are supported in {.arg arms}.")
    }
    arm_values <- factor(arm_values, arm_levels)
    arm_values <- droplevels(arm_values)

    # we split ref_ranges by `chr` and `arm` and then combine ranges in
    # each groups if there aren't any intervals.
    split_data <- data.table::data.table(
        chr = as.character(GenomeInfoDb::seqnames(ref_ranges)),
        arm = arm_values
    )

    # we create factor levels to determine proper order, we order chr_arm pairs
    # by chromosome firstly and then by arm.
    # for chr, we order it by integer portion and then by character
    chr_arm_pair <- data.table::copy(split_data)
    chr_arm_pair <- unique(chr_arm_pair)
    # nolint start
    chr_arm_pair[, seq_chr := sub("^chr", "", chr, perl = TRUE)]
    suppressWarnings(chr_arm_pair[, seq_int := as.integer(seq_chr)])
    chr_arm_pair[, chr_arm_order := order(seq_int, seq_chr, arm)]
    chr_arm_levels <- chr_arm_pair[, paste0(chr, arm)[chr_arm_order]]
    # nolint end

    split_factor <- factor(
        paste0(split_data[["chr"]], split_data[["arm"]]),
        chr_arm_levels
    )
    ref_ranges <- GenomicRanges::split(ref_ranges, split_factor, drop = TRUE)
    arm_gr_ranges <- S4Vectors::endoapply(ref_ranges, function(chr_cytoband) {
        gr <- GenomicRanges::reduce(chr_cytoband)
        S4Vectors::mcols(gr)[[arm_field]] <- unique(
            S4Vectors::mcols(chr_cytoband)[[arm_field]]
        )
        gr
    })
    if (any(lengths(arm_gr_ranges) > 1L)) {
        cli::cli_warn(
            "Cannot combine all ranges into one arm-level ranges",
            "i" = "Please check if {.arg ref_ranges} has intervals"
        )
    }
    arm_gr_ranges <- unlist(arm_gr_ranges, use.names = TRUE)
    arm_gr_ranges[S4Vectors::mcols(arm_gr_ranges)[[arm_field]] %in% arms]
}

#' Get UCSC cytoband data
#' @param x One of "hg19" or "hg38". "hg38" is derived from AnnotationHub by
#'   record id "AH53178" and "hg19" is by record id "AH53177". Default: "hg38".
#' @param add_arm A scalar logical value indicates whether add a column named
#'  "arm" defining the chromosome-arm for each items.
#' @return A [GenomicRanges][GenomicRanges::GRanges-class] object containing
#'   cytoband informations.
#' @export
get_cytoband <- function(x = "hg38", add_arm = TRUE) {
    out <- switch(x,
        hg19 = run_arm_cnv_ref_cytoband_hg19, # nolint
        hg38 = run_arm_cnv_ref_cytoband_hg38 # nolint
    )
    if (add_arm) {
        S4Vectors::mcols(out)$arm <- factor(
            data.table::fifelse(
                S4Vectors::mcols(out)$gieStain == "acen", "acen",
                sub("^([pq])[0-9.]+", "\\1", S4Vectors::mcols(out)$name)
            ),
            levels = c("p", "acen", "q", "")
        )
    }
    out
}

seg_to_arm <- function(seg_cnv, arm_cytoband, arm_field, group_fields = NULL, other_fields = character()) {
    # find ouverlap index --------------------------------
    overlap_hits <- silent_expr(
        GenomicRanges::findOverlaps(
            seg_cnv, arm_cytoband,
            type = "any"
        ),
        "Each of the 2 combined objects has sequence levels not in the other:",
        fixed = TRUE
    )
    if (length(overlap_hits) == 0L) {
        cli::cli_abort("No match between {.arg seg_cnv} and {.arg arm_cytoband}")
    }
    seg_hits <- S4Vectors::queryHits(overlap_hits)
    cytoband_hits <- S4Vectors::subjectHits(overlap_hits)

    # extract intersected ranges ------------------------
    intersect_region <- silent_expr(
        GenomicRanges::pintersect(
            seg_cnv[seg_hits], arm_cytoband[cytoband_hits],
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
        other_data <- data.table::as.data.table(seg_cnv[seg_hits])
        other_data <- other_data[, .SD, .SDcols = other_fields]
    }
    out <- cbind(out, other_data)

    # ensure every reference arm has a value
    ref_cytoband_dt <- unique(data.table::as.data.table(arm_cytoband))
    data.table::setnames(
        ref_cytoband_dt,
        c("seqnames", arm_field, "width"), c("chr", "arm", "arm_width")
    )
    ref_cytoband_dt <- ref_cytoband_dt[, c("chr", "arm", "arm_width")]
    out <- out[,
        .SD[ref_cytoband_dt, on = c("chr", "arm"), allow.cartesian = FALSE],
        by = group_fields
    ]
    # finally, move arm_width next to arm column
    data.table::setcolorder(out, "arm_width", after = "arm")
}

map_seqnames <- function(x, style, arg = rlang::caller_arg(x)) {
    if (all(GenomeInfoDb::seqlevelsStyle(x) != style)) {
        cli::cli_inform("Mapping seqnames of {.arg {arg}} to {style} style")
        GenomeInfoDb::seqlevelsStyle(x) <- style
    }
    x
}

assert_range_unique <- function(gr, group = NULL, arg_group = rlang::caller_arg(group), call = parent.frame()) {
    if (is.null(group)) {
        if (any_overlap(gr)) {
            cli::cli_abort(c(
                "Find overlapped ranges",
                i = "try to set {.arg {arg_group}}"
            ), call = call)
        }
    } else {
        gr_list <- split(gr, S4Vectors::mcols(gr)[[group]], drop = TRUE)
        failed_samples <- names(gr_list)[
            vapply(gr_list, any_overlap, logical(1L))
        ]
        if (length(failed_samples)) {
            cli::cli_abort(
                "Find overlapped ranges in group{?s}: {.val {failed_samples}}",
                call = call
            )
        }
    }
}

any_overlap <- function(gr) {
    hits <- GenomicRanges::findOverlaps(gr, gr)
    any(S4Vectors::queryHits(hits) != S4Vectors::subjectHits(hits))
}

gr_extend <- function(x, extension = 1L, use.names = TRUE) {
    GenomicRanges::update_ranges(x,
        start = GenomicRanges::start(x) - extension,
        end = GenomicRanges::end(x) + extension,
        use.names = use.names
    )
}

prepare_granges <- function(data, chr_field, start_field, end_field, other_fields = NULL, keep.extra.columns = TRUE, ignore.strand = TRUE, fun = NULL, call = parent.frame()) {
    assert_pkg("S4Vectors", fun = fun, frame = call)
    assert_pkg("GenomicRanges", fun = fun, frame = call)
    assert_class(chr_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        arg = rlang::caller_arg(chr_field),
        cross_msg = NULL, call = call
    )
    assert_class(start_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        arg = rlang::caller_arg(start_field),
        cross_msg = NULL, call = call
    )
    assert_class(end_field, rlang::is_scalar_character,
        "scalar {.cls character}",
        arg = rlang::caller_arg(end_field),
        cross_msg = NULL, call = call
    )
    assert_df_with_columns(data,
        c(other_fields, chr_field, start_field, end_field),
        arg = rlang::caller_arg(data), call = call
    )
    GenomicRanges::makeGRangesFromDataFrame(
        df = as.data.frame(data,
            check.names = FALSE,
            make.names = FALSE,
            stringsAsFactors = FALSE
        ),
        keep.extra.columns = keep.extra.columns,
        seqnames.field = chr_field,
        start.field = start_field,
        end.field = end_field,
        ignore.strand = ignore.strand
    )
}
