make_snv_nucleotide <- function(extension = 1L) {
    pairs_list <- c(
        rep_len(list(Biostrings::DNA_BASES), extension), 
        list(c("T", "C")),
        rep_len(list(Biostrings::DNA_BASES), extension)
    )
    n_pairs <- expand.grid(pairs_list)
    n_nucleotide <- .mapply(paste, n_pairs, list(sep = ""))
    unlist(n_nucleotide, recursive = FALSE, use.names = FALSE)
}

standardize_substitution <- function(x) {
    unname(substitution_pairs[x])
}

# values should be the standardized substitution
# names are substitution can be found out
#' @keywords internal
#' @noRd 
substitution_pairs <- structure(
    c(
        "T>C", "T>C", "C>T", "C>T", "T>A", "T>A", "T>G", "T>G",
        "C>A", "C>A", "C>G", "C>G"
    ),
    names = c(
        "A>G", "T>C", "C>T", "G>A", "A>T", "T>A", "A>C", "T>G",
        "C>A", "G>T", "C>G", "G>C"
    )
)

define_snv_motif <- function(motif, ref, alt, motif_extension = 1L) {
    substitution <- paste0(ref, ">", alt)

    substitution_motif <- combine_snv_motif(
        substitution, motif, motif_extension
    )

    # need to reverse-complement triplet for mutated purines (not just the
    # middle base) which ones need to be reverse-complemented
    substitution_type <- standardize_substitution(substitution)
    substitution_type_motif <- substitution_motif
    need_complement <- which(ref %chin% c("G", "A"))
    complemented_motif <- Biostrings::reverseComplement(motif)
    substitution_type_motif[need_complement] <- combine_snv_motif(
        substitution_type[need_complement],
        complemented_motif[need_complement],
        motif_extension
    )
    data.table::data.table(
        Substitution = factor(substitution,
            levels = names(substitution_pairs)
        ),
        SubstitutionMotif = factor(substitution_motif,
            levels = snv_motif
        ),
        SubstitutionType = factor(substitution_type,
            levels = unique(substitution_pairs)
        ),
        SubstitutionTypeMotif = factor(substitution_type_motif,
            levels = snv_type_motif
        )
    )
}

snv_motif <- c(
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
snv_type_motif <- c(
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

combine_snv_motif <- function(substitution, motif, motif_extension = 1L) {
    paste0(
        Biostrings::substring(motif, 1L, motif_extension),
        "[", substitution, "]",
        Biostrings::substring(
            motif, motif_extension + 2L,
            motif_extension * 2L + 1L
        )
    )
}
