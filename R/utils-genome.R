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
                "i" = "try to set {.arg arm_field} manually"
            ))
        }
        arm_field <- arm_field[[1L]]
    } else if (!rlang::is_scalar_character(arm_field)) {
        cli::cli_abort("{.arg arm_field} should be a scalar string")
    }
    arm_values <- S4Vectors::mcols(ref_ranges)[[arm_field]]
    arm_levels <- c("p", "acen", "q", "")
    assert_in(arm_values, arm_levels,
        arg_x = sprintf("{.arg %s}", arm_field)
    )
    assert_in(arms, arm_levels)
    arm_values <- factor(arm_values, arm_levels)
    arm_values <- droplevels(arm_values)

    # we split ref_ranges by `chr` and `arm` and then combine ranges in
    # each groups if there aren't any intervals.
    chr_values <- as.character(GenomeInfoDb::seqnames(ref_ranges))
    split_data <- paste0(chr_values, arm_values)
    split_data <- factor(
        split_data,
        unique(split_data[
            order(GenomeInfoDb::rankSeqlevels(chr_values), arm_values)
        ])
    )

    ref_ranges <- GenomicRanges::split(ref_ranges, split_data, drop = TRUE)
    out <- S4Vectors::endoapply(ref_ranges, function(chr_cytoband) {
        gr <- GenomicRanges::reduce(chr_cytoband)
        S4Vectors::mcols(gr)[[arm_field]] <- unique(
            S4Vectors::mcols(chr_cytoband)[[arm_field]]
        )
        gr
    })
    if (any(lengths(out) > 1L)) {
        cli::cli_warn(
            "Cannot combine all ranges into one arm-level ranges",
            "i" = "Please check if {.arg ref_ranges} has intervals"
        )
    }
    out <- unlist(out, use.names = TRUE)
    out[S4Vectors::mcols(out)[[arm_field]] %in% arms]
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
                sub("^([pq])[0-9.]+", "\\1",
                    S4Vectors::mcols(out)$name,
                    perl = TRUE
                )
            ),
            levels = c("p", "acen", "q", "")
        )
    }
    out
}

map_seqnames <- function(x, style, arg = rlang::caller_arg(x), style_arg = NULL) {
    if (all(GenomeInfoDb::seqlevelsStyle(x) != style)) {
        msg <- "Mapping seqnames of {.arg {arg}} to {style} style"
        if (!is.null(style_arg)) {
            msg <- paste(msg, "of {.arg {style_arg}}")
        }
        cli::cli_inform(msg)
        GenomeInfoDb::seqlevelsStyle(x) <- style
    }
    x
}


#' Get available/installed genomes
#' @param ref_genome BSgenome object or name of the installed BSgenome package.
#'  Default: "hg19". Details see "genome" argument of
#'  [getBSgenome][BSgenome::getBSgenome].
#' @keywords internal
get_genome <- function(ref_genome) {
    assert_pkg("BSgenome")
    ref_genome <- ref_genome %||% "hg19"
    BSgenome::getBSgenome(ref_genome)
}

assert_range_unique <- function(gr, group = NULL, arg_group = rlang::caller_arg(group), call = parent.frame()) {
    if (is.null(group)) {
        if (granges_any_overlap(gr)) {
            cli::cli_abort(c(
                "Find overlapped ranges",
                i = "try to set {.arg {arg_group}}"
            ), call = call)
        }
    } else {
        gr_list <- split(gr, S4Vectors::mcols(gr)[[group]], drop = TRUE)
        failed_samples <- names(gr_list)[
            vapply(gr_list, granges_any_overlap, logical(1L))
        ]
        if (length(failed_samples)) {
            cli::cli_abort(
                "Find overlapped ranges in group{?s}: {.val {failed_samples}}",
                call = call
            )
        }
    }
}

granges_any_overlap <- function(gr) {
    hits <- GenomicRanges::findOverlaps(gr, gr)
    any(S4Vectors::queryHits(hits) != S4Vectors::subjectHits(hits))
}

granges_extend <- function(x, extension = 1L, use.names = TRUE) {
    GenomicRanges::update_ranges(x,
        start = GenomicRanges::start(x) - extension,
        end = GenomicRanges::end(x) + extension,
        use.names = use.names
    )
}
