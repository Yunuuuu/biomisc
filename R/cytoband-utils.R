#' Get arm-level ranges
#' @param ref_ranges a [`GenomicRanges`][GenomicRanges::GRanges-class] object to
#' combine into arm-level ranges
#' @param arm_col a scalar string indicates the column containing the chromosome
#' arm
#' @return a [`GenomicRanges`][GenomicRanges::GRanges-class] object containing
#' arm-level ranges.
#' @export 
arm_ranges <- function(ref_ranges, arm_col = NULL) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
        cli::cli_abort("GenomicRanges must be installed to use this function.",
            call. = FALSE
        )
    }
    if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
        cli::cli_abort("GenomeInfoDb must be installed to use this function.")
    }
    if (!requireNamespace("S4Vectors", quietly = TRUE)) {
        cli::cli_abort("S4Vectors must be installed to use this function.",
            call. = FALSE
        )
    }

    if (!inherits(ref_ranges, "GenomicRanges")) {
        cli::cli_abort("{.arg ref_ranges} should be a {.cls GenomicRanges} object")
    }

    if (!identical(GenomeInfoDb::seqlevelsStyle(ref_ranges), "UCSC")) {
        cli::cli_inform("try to map seqnames of ref_ranges to UCSC style")
        GenomeInfoDb::seqlevels(ref_ranges) <- "UCSC"
    }
    if (is.null(arm_col)) {
        arm_col <- grep("^arm", names(S4Vectors::mcols(ref_ranges)),
            ignore.case = TRUE,
            value = TRUE
        )
        arm_col <- arm_col[[1L]]
        if (!length(arm_col)) {
            cli::cli_abort(c(
                "Cannot determine right {.arg arm_col}",
                "i" = "try to set {.arg arm_col}"
            ))
        }
    } else if (!rlang::is_scalar_character(arm_col)) {
        cli::cli_abort("{.arg arm_col} should be a scalar string")
    }

    # we split ref_ranges by `chr` and `arm` and then combine ranges in
    # each groups if there aren't any intervals. 
    # we create factor levels to determine proper order, we order chr_arm pairs
    # by chromosome firstly and then by arm. 
    # for chr, we order it by integer portion and then by character
    chr_arm_dt <- data.table::data.table(
        chr = as.character(GenomeInfoDb::seqnames(ref_ranges)),
        arm = S4Vectors::mcols(ref_ranges)[[arm_col]]
    )
    chr_arm_pair <- unique(chr_arm_dt)
    chr_arm_pair[, seq_chr := sub("chr", "", chr)]
    suppressWarnings(
        chr_arm_pair[, seq_int := as.integer(seq_chr)]
    )
    chr_arm_pair[, chr_arm_order := order(seq_int, seq_chr, arm)]
    chr_arm_levels <- paste0(chr_arm_pair[["chr"]], chr_arm_pair[["arm"]])
    chr_arm_levels <- chr_arm_levels[chr_arm_pair[["chr_arm_order"]]]
    chr_arm_group <- factor(
        paste0(chr_arm_dt[["chr"]], chr_arm_dt[["arm"]]),
        chr_arm_levels
    )
    ref_ranges <- GenomicRanges::split(
        ref_ranges, chr_arm_group,
        drop = TRUE
    )
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
