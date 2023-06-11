# for mutation (including SNV, DBS, MBS)
# for substitution (sub), we indicates with string ">"
combine_mut_context <- function(context, sub, extension = 1L) {
    paste0(
        Biostrings::substring(context, 1L, extension),
        "[", sub, "]",
        Biostrings::substring(
            context, Biostrings::nchar(context) - extension,
            Biostrings::nchar(context)
        )
    )
}

standard_snv_context <- function(extension = 1L) {
    pairs_list <- c(
        rep_len(list(Biostrings::DNA_BASES), extension),
        list(c("T", "C")),
        rep_len(list(Biostrings::DNA_BASES), extension)
    )
    n_pairs <- expand.grid(pairs_list)
    n_nucleotide <- .mapply(paste, n_pairs, list(sep = ""))
    unlist(n_nucleotide, recursive = FALSE, use.names = FALSE)
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

standardize_snv_sub <- function(x) {
    unname(standard_snv_sub[x])
}

sub_context <- function(chr, pos, ref, alt, strand = "+", ref_genome = NULL, extension = 1L, bg_extension = NULL) {
    # get ref_genome
    ref_genome <- ref_genome %||% "hg19"
    ref_genome <- BSgenome::getBSgenome(ref_genome)

    # create mutation GenomicRanges
    mut_gr <- GenomicRanges::makeGRangesFromDataFrame(
        data.frame(
            chr = chr, start = pos,
            end = pos + nchar(ref) - 1L,
            ref = ref, alt = alt, strand = strand
        ),
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end",
        strand.field = "strand",
        starts.in.df.are.0based = FALSE
    )

    # match seqlevelstyle between mut_gr and ref_genome
    ref_genome_style <- GenomeInfoDb::seqlevelsStyle(ref_genome)
    mut_gr <- map_seqnames(mut_gr, ref_genome_style, "chr")
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

    # Get context
    context <- suppressWarnings(gr_extend(mut_gr, extension = extension))
    good_ranges <- GenomicRanges::trim(context) == context
    mut_gr <- mut_gr[good_ranges]
    context <- BSgenome::getSeq(
        x = ref_genome, context[good_ranges],
        as.character = FALSE
    )

    # check mutation reference allele match the second base in context
    matched_ranges <- mut_gr$ref == as.character(
        Biostrings::subseq(context, extension + 1L, extension + 1L)
    )
    if (any(!matched_ranges)) {
        cli::cli_warn("Removing {.val {sum(!matched_ranges)}} ({format_percent(mean(!matched_ranges))}) mutations whose reference allele cannot match the {.arg ref_genome}")
    }
    mut_gr <- mut_gr[matched_ranges]
    context <- context[matched_ranges]

    if (is.null(bg_extension)) {
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

    sub <- paste0(mut_gr$ref, ">", mut_gr$alt)
    sub_type <- standardize_snv_sub(sub)
    # need to reverse-complement triplet for mutated purines (not just the
    # middle base) which ones need to be reverse-complemented
    need_complement <- sub_type != sub
    sub_type_context <- context
    sub_type_context[need_complement] <- Biostrings::reverseComplement(
        sub_type_context[need_complement]
    )
    list(
        mut = mut_gr, sub = sub,
        sub_type = sub_type,
        context = context,
        sub_type_context = sub_type_context,
        bg_context = bg
    )
}

define_snv_sub_tri_context <- function(motif, ref, alt) {
    sub <- paste0(ref, ">", alt)

    sub_motif <- combine_mut_context(
        motif, sub, 1L
    )

    # need to reverse-complement triplet for mutated purines (not just the
    # middle base) which ones need to be reverse-complemented
    sub_type <- standardize_snv_sub(sub)
    sub_type_motif <- sub_motif
    need_complement <- which(ref %chin% c("G", "A"))
    sub_type_motif[need_complement] <- combine_mut_context(
        Biostrings::reverseComplement(motif)[need_complement],
        sub_type[need_complement], 1L
    )
    data.table::data.table(
        Substitution = factor(sub,
            levels = names(standard_snv_sub)
        ),
        SubstitutionMotif = factor(sub_motif,
            levels = snv_sub_tri_context
        ),
        SubstitutionType = factor(sub_type,
            levels = unique(standard_snv_sub)
        ),
        SubstitutionTypeMotif = factor(sub_type_motif,
            levels = standard_snv_sub_tri_context
        )
    )
}

snv_sub_tri_context <- c(
    "C[G>A]A", "T[C>G]C", "A[G>A]C", "T[C>T]A", "T[G>C]A", "C[G>A]C",
    "A[A>G]C", "T[G>A]G", "G[G>A]A", "G[C>T]G", "T[G>A]A", "T[C>A]T",
    "C[C>T]A", "G[G>C]A", "T[C>A]C", "T[C>T]G", "A[G>C]A", "G[C>T]T",
    "G[A>G]G", "T[C>T]T", "C[C>T]C", "C[G>A]G", "T[T>A]C", "C[G>A]T",
    "C[C>T]G", "A[G>A]A", "T[G>T]G", "A[G>T]A", "C[C>G]T", "A[T>C]A",
    "G[G>T]A", "T[C>T]C", "T[C>A]A", "T[C>G]A", "A[T>G]C", "G[G>A]C",
    "A[G>T]T", "A[G>A]G", "T[C>G]T", "C[C>A]A", "C[T>G]A", "C[A>G]C",
    "C[C>A]C", "T[C>G]G", "A[C>G]G", "C[C>G]C", "T[G>T]T", "T[G>T]C",
    "G[G>T]C", "A[T>C]G", "C[T>G]T", "A[C>T]G", "G[T>C]A", "G[A>C]G",
    "A[G>C]T", "T[G>C]G", "A[C>G]A", "C[G>C]A", "T[G>C]T", "T[G>A]C",
    "A[T>C]C", "T[G>A]T", "A[T>G]T", "C[T>A]T", "A[A>G]T", "A[T>A]G",
    "G[G>T]G", "G[G>C]C", "G[C>A]T", "T[G>T]A", "A[C>A]A", "C[C>A]T",
    "C[T>C]C", "G[G>A]G", "G[A>C]A", "A[T>C]T", "T[A>T]T", "C[A>T]G",
    "A[C>T]C", "T[T>C]C", "G[C>T]C", "A[A>G]A", "A[C>A]T", "G[T>C]C",
    "A[G>T]G", "A[T>A]T", "G[G>C]T", "C[G>T]G", "A[C>G]C", "C[G>T]A",
    "A[C>T]A", "C[A>T]T", "C[G>T]T", "A[G>T]C", "G[C>T]A", "C[C>A]G",
    "G[C>A]A", "G[C>A]G", "G[T>G]C", "G[G>T]T", "G[A>G]C", "T[C>A]G",
    "C[C>T]T", "A[G>A]T", "A[A>T]C", "C[T>C]A", "C[C>G]A", "T[T>A]T",
    "C[T>C]G", "A[G>C]C", "G[C>G]C", "A[C>T]T", "T[G>C]C", "G[C>G]T",
    "G[C>A]C", "C[T>A]C", "C[A>G]T", "A[A>C]T", "G[G>C]G", "A[G>C]G",
    "T[T>C]T", "A[A>T]G", "C[G>C]C", "C[T>A]G", "T[T>A]A", "T[A>C]T",
    "A[C>A]C", "G[T>C]T", "A[T>A]C", "C[A>C]T", "T[A>C]C", "G[A>T]G",
    "T[A>C]A", "G[A>G]T", "G[A>T]T", "C[T>C]T", "G[T>A]T", "G[A>G]A",
    "T[T>G]A", "T[T>G]T", "T[A>G]T", "G[A>T]C", "A[A>C]A", "G[G>A]T",
    "C[A>C]G", "C[A>G]A", "G[T>G]A", "G[A>T]A", "T[A>C]G", "A[C>G]T",
    "G[T>A]A", "T[A>G]C", "C[G>C]G", "T[T>C]G", "T[A>G]G", "A[A>T]T",
    "T[T>G]G", "C[A>T]A", "C[A>C]C", "C[G>T]C", "G[T>A]C", "C[A>G]G",
    "C[A>C]A", "T[A>T]A", "G[T>C]G", "A[A>C]C", "C[A>T]C", "C[T>G]G",
    "A[A>G]G", "G[C>G]A", "A[A>C]G", "T[A>T]G", "C[G>C]T", "T[T>G]C",
    "A[T>G]G", "G[T>G]T", "G[T>A]G", "C[C>G]G", "A[C>A]G", "A[A>T]A",
    "A[T>G]A", "G[A>C]T", "G[T>G]G", "G[A>C]C", "A[T>A]A", "C[T>G]C",
    "T[A>G]A", "G[C>G]G", "T[T>C]A", "C[T>A]A", "T[T>A]G", "T[A>T]C"
)

# Possible substitution types after being referred to by the pyrimidine of the mutated Watson-Crick base pair
standard_snv_sub_tri_context <- c(
    "A[C>A]A", "A[C>A]C", "A[C>A]G", "A[C>A]T", "C[C>A]A", "C[C>A]C",
    "C[C>A]G", "C[C>A]T", "G[C>A]A", "G[C>A]C", "G[C>A]G", "G[C>A]T",
    "T[C>A]A", "T[C>A]C", "T[C>A]G", "T[C>A]T", "A[C>G]A", "A[C>G]C",
    "A[C>G]G", "A[C>G]T", "C[C>G]A", "C[C>G]C", "C[C>G]G", "C[C>G]T",
    "G[C>G]A", "G[C>G]C", "G[C>G]G", "G[C>G]T", "T[C>G]A", "T[C>G]C",
    "T[C>G]G", "T[C>G]T", "A[C>T]A", "A[C>T]C", "A[C>T]G", "A[C>T]T",
    "C[C>T]A", "C[C>T]C", "C[C>T]G", "C[C>T]T", "G[C>T]A", "G[C>T]C",
    "G[C>T]G", "G[C>T]T", "T[C>T]A", "T[C>T]C", "T[C>T]G", "T[C>T]T",
    "A[T>A]A", "A[T>A]C", "A[T>A]G", "A[T>A]T", "C[T>A]A", "C[T>A]C",
    "C[T>A]G", "C[T>A]T", "G[T>A]A", "G[T>A]C", "G[T>A]G", "G[T>A]T",
    "T[T>A]A", "T[T>A]C", "T[T>A]G", "T[T>A]T", "A[T>C]A", "A[T>C]C",
    "A[T>C]G", "A[T>C]T", "C[T>C]A", "C[T>C]C", "C[T>C]G", "C[T>C]T",
    "G[T>C]A", "G[T>C]C", "G[T>C]G", "G[T>C]T", "T[T>C]A", "T[T>C]C",
    "T[T>C]G", "T[T>C]T", "A[T>G]A", "A[T>G]C", "A[T>G]G", "A[T>G]T",
    "C[T>G]A", "C[T>G]C", "C[T>G]G", "C[T>G]T", "G[T>G]A", "G[T>G]C",
    "G[T>G]G", "G[T>G]T", "T[T>G]A", "T[T>G]C", "T[T>G]G", "T[T>G]T"
)
