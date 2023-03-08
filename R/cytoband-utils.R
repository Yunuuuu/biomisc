#' Get arm-level ranges
#' @param ref_ranges A [`GenomicRanges`][GenomicRanges::GRanges-class] object to
#'   combine into arm-level ranges.
#' @param arm_col A scalar string indicates the column containing the chromosome
#'   arm. If `NULL`, the internal will search the first column starts with "arm"
#'   (ignore letter case). Only values of "p", "q" and "acen" are supported.
#' @return A [GenomicRanges][GenomicRanges::GRanges-class] object containing
#' arm-level ranges.
#' @export
get_arm_ranges <- function(ref_ranges, arm_col = NULL) {
    assert_pkg("S4Vectors")
    assert_pkg("GenomicRanges")
    assert_pkg("GenomeInfoDb")
    assert_class(ref_ranges, "GenomicRanges")

    if (!identical(GenomeInfoDb::seqlevelsStyle(ref_ranges), "UCSC")) {
        cli::cli_inform("try to map seqnames of {.arg ref_ranges} to UCSC style")
        GenomeInfoDb::seqlevels(ref_ranges) <- "UCSC"
    }
    if (is.null(arm_col)) {
        arm_col <- grep("^arm", names(S4Vectors::mcols(ref_ranges)),
            ignore.case = TRUE, perl = TRUE,
            value = TRUE
        )
        if (!length(arm_col)) {
            cli::cli_abort(c(
                "Cannot determine right {.arg arm_col}",
                "i" = "try to set {.arg arm_col}"
            ))
        }
        arm_col <- arm_col[[1L]]
    } else if (!rlang::is_scalar_character(arm_col)) {
        cli::cli_abort("{.arg arm_col} should be a scalar string")
    }
    arm_values <- S4Vectors::mcols(ref_ranges)[[arm_col]]
    arm_levels <- c("p", "acen", "q")
    if (!any(arm_values %in% arm_levels)) {
        cli::cli_abort("Only values of {.val {arm_levels}} are supported in column specified in {.arg arm_col}.")
    }
    arm_values <- factor(arm_values, arm_levels)
    arm_values <- droplevels(arm_values)

    # we split ref_ranges by `chr` and `arm` and then combine ranges in
    # each groups if there aren't any intervals.
    split_data <- data.table::data.table(
        chr = as.character(GenomeInfoDb::seqnames(ref_ranges)),
        arm = factor(S4Vectors::mcols(ref_ranges)[[arm_col]], c(""))
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
    arm_ranges <- S4Vectors::endoapply(ref_ranges, function(chr_cytoband) {
        GenomicRanges::reduce(chr_cytoband)
    })
    if (any(lengths(arm_ranges) > 1L)) {
        cli::cli_warn(
            "Cannot combine all ranges into one arm-level ranges",
            "i" = "Please check if {.arg ref_ranges} has intervals"
        )
    }
    arm_ranges <- unlist(arm_ranges, use.names = TRUE)
    arm_ranges
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
            levels = c("p", "acen", "q")
        )
    }
    out
}
