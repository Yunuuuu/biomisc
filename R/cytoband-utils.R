#' Get arm-level ranges
#' @param ref_ranges a [`GenomicRanges`][GenomicRanges::GRanges-class] object to
#' combine into arm-level ranges
#' @param arm_col a scalar string indicates the column containing the chromosome
#' arm
#' @return a [`GenomicRanges`][GenomicRanges::GRanges-class] object containing
#' arm-level ranges.
#' @export
get_arm_ranges <- function(ref_ranges, arm_col = NULL) {
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

#' Prepare data for the input of PyClone or PyClone-vi
#' @param mut_data a data.frame of mutation data with the top four columns
#' containing "chromosome", "pos", "ref_counts" and "var_counts", other columns
#' will be omited. Column names don't matter.
#' @param cnv_data a data.frame of allelic CNV data with the top five columns
#' containing "chromosome", "start_pos", "end_pos", "major_cn", and "minor_cn",
#' other columns will be omited. Column names don't matter.
#' @param normal_cn The copy number of the locus in non-malignant cells. This
#' should generally be 2 except for sex chromosomes in males. 
#' @param pyclone_vi a logical value indicates whether prepare data for the
#' input of PyClone-vi. Details see <https://github.com/Roth-Lab/pyclone-vi>
#' @param sample_id Only used for
#' [PyClone-vi](https://github.com/Roth-Lab/pyclone-vi), the sample names.
#' @param purity Only used for
#' [PyClone-vi](https://github.com/Roth-Lab/pyclone-vi), the tumour content
#' (cellularity) of the sample. Default value is 1.0 if column is not present. 
#' @param error_rate Only used for
#' [PyClone-vi](https://github.com/Roth-Lab/pyclone-vi), sequencing error rate.
#' @references 
#'  - <https://github.com/Roth-Lab/pyclone>
#'  - <https://github.com/Roth-Lab/pyclone-vi>
#' @export 
prepare_pyclone <- function(mut_data, cnv_data, normal_cn = 2, pyclone_vi = FALSE, sample_id = NULL, purity = 1L, error_rate = NULL) {
    if (pyclone_vi) {
        if (!rlang::is_scalar_character(sample_id)) {
            cli::cli_abort(
                "{.arg sample_id} should be a scalar string"
            )
        }
        if (!(rlang::is_scalar_double(purity) || rlang::is_scalar_integer(purity))) {
            # nolint
            cli::cli_abort("{.arg purity} should be a scalar numeric")
        }
    }
    mut_data <- data.table::as.data.table(mut_data)[, 1:4]
    data.table::setnames(
        mut_data, c("chromosome", "pos", "ref_counts", "var_counts")
    )

    cnv_data <- data.table::as.data.table(cnv_data)[, 1:5]
    data.table::setnames(
        cnv_data,
        c("chromosome", "start_pos", "end_pos", "major_cn", "minor_cn")
    )

    pyclone_data <- cnv_data[
        mut_data,
        list(
            ref_counts = i.ref_counts, var_counts = i.var_counts,
            minor_cn = x.minor_cn, major_cn = x.major_cn,
            chromosome = i.chromosome, position = i.pos,
            start_pos = x.start_pos, end_pos = x.end_pos
        ),
        on = c("chromosome", "start_pos<=pos", "end_pos>=pos"),
        nomatch = NULL
    ][, .SD[rowSums(is.na(.SD)) == 0L]]
    test_pos <- pyclone_data[["pos"]] < pyclone_data[["start_pos"]] |
        pyclone_data[["pos"]] > pyclone_data[["end_pos"]]
    if (any(test_pos)) {
        cli::cli_abort("Something wrong when parsing CN of mutation")
    }
    pyclone_data[, mutation_id := paste(chromosome, position, sep = ":")]
    pyclone_data[, normal_cn := normal_cn]

    data.table::setcolorder(
        pyclone_data,
        c(
            "mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn",
            "major_cn"
        )
    )
    if (pyclone_vi) {
        pyclone_data[, sample_id := sample_id]
        pyclone_data[, tumour_content := purity]
        data.table::setnames(pyclone_data, "var_counts", "alt_counts")
        data.table::setcolorder(pyclone_data,
            "sample_id",
            after = "mutation_id"
        )
        data.table::setcolorder(pyclone_data,
            c("minor_cn", "normal_cn"),
            after = "major_cn"
        )
        data.table::setcolorder(pyclone_data,
            "tumour_content",
            after = "normal_cn"
        )
        if (!is.null(error_rate)) {
            pyclone_data[, error_rate := error_rate]
            data.table::setcolorder(pyclone_data,
                "error_rate",
                after = "tumour_content"
            )
        }
    }
    pyclone_data[major_cn > 0L]
}
