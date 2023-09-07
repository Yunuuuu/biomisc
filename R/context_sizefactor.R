#' Evaluate the sizefactor to transform context between different scale
#' @param context A character of nucleotide sequence, all must have the same
#'  number of characters.
#' @inheritParams get_genome
#' @param method One of "genome", "exome", "exome2genome", "genome2exome", If
#'  the input mut_data only contains the counts of the mutations observed in
#'  each context, then the mut_data must be normalized.
#' - If method is set to 'exome', the input data frame is normalized by number
#'   of times each trinucleotide context is observed in the exome.
#' - If method is set to 'genome', the input data frame is normalized by number
#'   of times each trinucleotide context is observed in the genome.
#' - If method is set to 'exome2genome', normalization is performed to reflect
#'   the absolute frequency of each trinucleotide context as it would occur
#'   across the whole genome. Thus the count data for each trinucleotide context
#'   is multiplied by a ratio of that trinucleotide's occurence in the genome to
#'   the trinucleotide's occurence in the exome.
#'
#'   Details see <https://github.com/raerose01/deconstructSigs>
#' @param exome_ranges A [GenomicRanges][GenomicRanges::GRanges] Object define
#'  exome ranges.
#' @param contigs An atomic vector specifying the chromosome to define the  size
#' factor to normalize the context matrix.
#' @param strand A character vector indicates the strand in which to calculate
#' the context frequency. Both positive and negtive strands.
#' @param pruning_mode When some of the seqlevels to drop from `exome_ranges`
#'  are in use (i.e. have ranges on them), the ranges on these sequences need to
#'  be removed before the seqlevels can be dropped. We call this pruning. The
#'  pruning_mode argument controls how to prune x. See
#'  [seqinfo][GenomeInfoDb::seqinfo] pruning.mode argument.
#' @return A numeric size factor for each mut_context
#' @export
context_sizefactor <- function(context, method = c("genome", "exome", "exome2genome", "genome2exome"), ref_genome = NULL, exome_ranges = NULL, contigs = NULL, strand = c("+", "-"), pruning_mode = "error") {
    assert_(
        context, function(x) {
            n <- nchar(x)
            all(n == n[1L])
        },
        what = "an atomic character with equal strings"
    )
    assert_inclusive(strand, c("+", "-"))
    ref_genome <- get_genome(ref_genome)
    calculate_context_sizefactor(
        context, ref_genome,
        exome_ranges = exome_ranges,
        method = method,
        contigs = contigs,
        strand = unique(strand),
        pruning_mode = pruning_mode
    )
}

calculate_context_sizefactor <- function(context, ref_genome = NULL, exome_ranges = NULL, method = NULL, contigs = NULL, strand = c("+", "-"), pruning_mode = "error") {
    method <- match.arg(
        method,
        c("genome", "exome", "exome2genome", "genome2exome")
    )
    context_width <- nchar(context)[1L]
    if (grepl("genome", method, fixed = TRUE)) {
        seqstyle <- GenomeInfoDb::seqlevelsStyle(ref_genome)
        seqstyle_arg <- "ref_genome"
    } else {
        seqstyle <- GenomeInfoDb::seqlevelsStyle(exome_ranges)
        seqstyle_arg <- "exome_ranges"
    }
    if (length(seqstyle) > 1L) {
        cli::cli_abort(" the styles of the seqlevels in {.arg {seqstyle_arg}} is not unique")
    }
    if (!is.null(contigs)) {
        contigs <- as.character(contigs)
        contigs <- map_seqnames(contigs, seqstyle, style_arg = seqstyle_arg)
    }
    if (grepl("genome", method, fixed = TRUE)) {
        if (is.null(contigs)) {
            contigs <- GenomeInfoDb::seqnames(ref_genome)
        }
        counts_wgs <- lapply(strand, function(x) {
            counts <- Biostrings::oligonucleotideFrequency(
                BSgenome::getSeq(ref_genome, contigs, strand = x),
                width = context_width
            )
            colSums(counts)
        })
        wgs_contexts <- Reduce(intersect, lapply(counts_wgs, names))
        wgs_missing_contexts <- setdiff(context, wgs_contexts) # nolint
        if (length(wgs_missing_contexts) > 0L) {
            cli::cli_abort(c(
                "Cannot find all {.arg context} in {.arg ref_genome}",
                x = "missing context{?s}: {.val {wgs_missing_contexts}}"
            ))
        }
        counts_wgs <- Reduce(`+`, lapply(counts_wgs, `[`, context))
    }
    if (grepl("exome", method, fixed = TRUE)) {
        exome_ranges <- map_seqnames(exome_ranges,
            seqstyle,
            style_arg = seqstyle_arg
        )
        exome_ranges <- GenomicRanges::reduce(exome_ranges)
        if (!is.null(contigs)) {
            exome_ranges <- GenomeInfoDb::keepSeqlevels(
                exome_ranges, contigs,
                pruning.mode = pruning_mode
            )
        }
        counts_wes <- lapply(strand, function(x) {
            GenomicRanges::strand(exome_ranges) <- x
            counts <- Biostrings::oligonucleotideFrequency(
                BSgenome::getSeq(ref_genome, exome_ranges),
                width = context_width
            )
            colSums(counts)
        })
        wes_contexts <- Reduce(intersect, lapply(counts_wes, names))
        wes_missing_contexts <- setdiff(context, wes_contexts) # nolint
        if (length(wes_missing_contexts) > 0L) {
            cli::cli_abort(c(
                "Cannot find all {.arg context} in {.arg exome_ranges}",
                x = "missing context{?s}: {.val {wes_missing_contexts}}"
            ))
        }
        counts_wes <- Reduce(`+`, lapply(counts_wes, `[`, context))
    }
    switch(method,
        # return mut counts divided by number of times that trinucleotide context is observed in the genome
        genome = 1L / counts_wgs,

        # return mut counts divided by number of times that trinucleotide
        # context is observed in the exome
        exome = 1L / counts_wes,

        # use the ratio of WGS/WES to normalize the input data
        exome2genome = counts_wgs / counts_wes,

        # use the ratio of WES/WGS to normalize the input data
        genome2exome = counts_wes / counts_wgs
    )
}
