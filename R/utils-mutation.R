# for mutation (including SNV, DBS, MBS)
# for substitution (sub), we indicates with string ">"
combine_sub_context <- function(context, sub, extension = 1L) {
    paste0(
        Biostrings::substring(context, 1L, extension),
        "[", sub, "]",
        Biostrings::substring(
            context, Biostrings::nchar(context) - extension + 1L,
            Biostrings::nchar(context)
        )
    )
}

# values should be the standardized substitution
# names are substitution can be found out
#' @keywords internal
#' @noRd
standard_snv_sub <- structure(
    c(
        "T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G",
        "C>A", "C>A", "C>G", "C>G"
    ),
    names = c(
        "A>G", "T>C", "C>T", "G>A", "A>T", "T>A", "A>C", "T>G",
        "C>A", "G>T", "C>G", "G>C"
    )
)

enumerate_mut_context <- function(mut_pattern, extension = 1L) {
    pairs_list <- c(
        rep_len(list(Biostrings::DNA_BASES), extension),
        list(mut_pattern),
        rep_len(list(Biostrings::DNA_BASES), extension)
    )
    n_pairs <- expand.grid(pairs_list)
    n_nucleotide <- .mapply(paste, n_pairs, list(sep = ""))
    unlist(n_nucleotide, recursive = FALSE, use.names = FALSE)
}

enumerate_standard_snv_context <- function(extension = 1L) {
    enumerate_mut_context(c("T", "C"), extension = extension)
}

enumerate_standard_snv_sub_context <- function(extension = 1L) {
    enumerate_mut_context(
        paste0("[", unique(standard_snv_sub), "]"),
        extension = extension
    )
}

enumerate_snv_sub_context <- function(extension = 1L) {
    enumerate_mut_context(
        paste0("[", unique(names(standard_snv_sub)), "]"),
        extension = extension
    )
}

sub_context_to_mut_context <- function(sub_context) {
    gsub(
        sprintf(
            "[\\[\\]]|\\>[%s]+",
            paste0(Biostrings::DNA_BASES, collapse = "")
        ), "",
        sub_context,
        perl = TRUE
    )
}

standardize_snv_sub <- function(x) {
    unname(standard_snv_sub[x])
}

snv_sub_context <- function(mut_data, chr_field = "chr", mut_pos = "pos", ref_field = "ref", alt_field = "alt", strand = "+", ref_genome, extension = 1L, bg_extension = NULL) {
    # create mutation GenomicRanges
    mut_gr <- GenomicRanges::makeGRangesFromDataFrame(
        as.data.frame(mut_data, make.names = FALSE, check.names = FALSE),
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = chr_field,
        start.field = mut_pos,
        end.field = mut_pos,
        starts.in.df.are.0based = FALSE
    )
    mut_gr <- mut_gr[
        nchar(S4Vectors::mcols(mut_gr)[[ref_field]]) ==
            nchar(S4Vectors::mcols(mut_gr)[[alt_field]]) &
            nchar(S4Vectors::mcols(mut_gr)[[ref_field]]) == 1L
    ]
    GenomicRanges::strand(mut_gr) <- strand

    # match seqlevelstyle between mut_gr and ref_genome
    ref_genome_style <- GenomeInfoDb::seqlevelsStyle(ref_genome)
    mut_gr <- map_seqnames(mut_gr, ref_genome_style, "mut_data")
    GenomeInfoDb::seqlevels(mut_gr) <- GenomeInfoDb::seqlevels(ref_genome)
    suppressWarnings(
        GenomeInfoDb::seqinfo(mut_gr) <- GenomeInfoDb::seqinfo(ref_genome)
    )

    # remove out-of ranges mutations
    out_ranges <- GenomicRanges::trim(mut_gr) != mut_gr
    if (any(out_ranges)) {
        cli::cli_warn("Removing {.val {sum(out_ranges)}} ({format_percent(mean(out_ranges))}) out-of-bound mutations compared with {.arg ref_genome}")
    }
    mut_gr <- mut_gr[!out_ranges]

    # check mutation reference allele match the second base in context
    matched_ranges <- S4Vectors::mcols(mut_gr)[[ref_field]] == as.character(
        Biostrings::getSeq(ref_genome, mut_gr)
    )
    if (any(!matched_ranges)) {
        cli::cli_warn("Removing {.val {sum(!matched_ranges)}} ({format_percent(mean(!matched_ranges))}) mutations whose reference allele cannot match the {.arg ref_genome}")
    }
    mut_gr <- mut_gr[matched_ranges]

    # Get context
    context <- suppressWarnings(gr_extend(mut_gr, extension = extension))
    good_ranges <- GenomicRanges::trim(context) == context
    mut_gr <- mut_gr[good_ranges]
    context <- BSgenome::getSeq(x = ref_genome, context[good_ranges])
    if (!is.null(bg_extension)) {
        bg <- suppressWarnings(gr_extend(mut_gr, extension = bg_extension))
        # omit out-of-bound ranges
        trimed_bg <- GenomicRanges::trim(bg)
        good_ranges <- trimed_bg == bg
        bg <- BSgenome::getSeq(
            x = ref_genome, trimed_bg, as.character = FALSE
        )
    } else {
        bg <- NULL
    }

    sub <- paste0(
        S4Vectors::mcols(mut_gr)[[ref_field]], ">",
        S4Vectors::mcols(mut_gr)[[alt_field]]
    )
    standard_snv_sub <- standardize_snv_sub(sub)
    # need to reverse-complement triplet for mutated purines (not just the
    # middle base) which ones need to be reverse-complemented
    need_complement <- which(standard_snv_sub != sub)
    standard_snv_context <- context
    standard_snv_context[need_complement] <- Biostrings::reverseComplement(
        standard_snv_context[need_complement]
    )
    list(
        mut = mut_gr,
        sub = sub,
        context = context,
        standard_snv_sub = standard_snv_sub,
        standard_snv_context = standard_snv_context,
        bg_context = bg
    )
}

define_snv_sub_tri_context <- function(motif, ref, alt) {
    sub <- paste0(ref, ">", alt)

    sub_motif <- combine_sub_context(motif, sub, 1L)

    # need to reverse-complement triplet for mutated purines (not just the
    # middle base) which ones need to be reverse-complemented
    sub_type <- standardize_snv_sub(sub)
    sub_type_motif <- sub_motif
    need_complement <- which(sub_type != sub)
    sub_type_motif[need_complement] <- combine_sub_context(
        Biostrings::reverseComplement(motif)[need_complement],
        sub_type[need_complement], 1L
    )
    data.table::data.table(
        Substitution = factor(sub,
            levels = names(standard_snv_sub)
        ),
        SubstitutionMotif = factor(sub_motif,
            levels = enumerate_snv_sub_context()
        ),
        SubstitutionType = factor(sub_type,
            levels = unique(standard_snv_sub)
        ),
        SubstitutionTypeMotif = factor(sub_type_motif,
            levels = enumerate_standard_snv_sub_context()
        )
    )
}
