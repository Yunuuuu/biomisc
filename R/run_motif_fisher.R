#' Extract single 5' and 3' bases flanking the mutated site for de-novo
#' signature analysis. Also estimates APOBEC enrichment scores.
#' @details Extracts immediate 5' and 3' bases flanking the mutated site and
#' classifies them into 96 substitution classes.
#' Requires BSgenome data packages for sequence extraction.
#'
#' APOBEC Enrichment: Enrichment score is calculated using the same method
#' described by Roberts et al.
#'
#'    E = (MUTATIONS(tCw) * CONTEXT(CorG)) / (MUTATIONS(CorG) * CONTEXT(tCw))
#'
#'  where, MUTATIONS(tCw) = number of mutations within T\[C>T\]W and T\[C>G\]W
#'  context. (W -> A or T)
#'
#'         MUTATIONS(CorG) = number of mutated C and G
#'
#' One-sided Fisher's Exact test is performed to determine the enrichment of
#' signature motif over background.
#' 
#' APOBEC_enrich mean the enrichment over random of APOBEC pattern mutations.
#' This is calculated as: {\[(TCW to TGW) + (TCW to TTW)\]/\[(C to G) + (C to
#' T)\]}/\[TCW/C\]. The enrichment value > 2, which implies that in such samples
#' at least 50% of APOBEC signature mutations have been in fact made by APOBEC
#' enzyme(s). APOBEC_MutLoad_MinEstimate is the minimum estimate of number of
#' APOBEC induced mutations in a sample. This estimate is calculated using the
#' formula: \[TCW to TGW + TCW to TTW\]×\[(APOBEC_enrich-1)/APOBEC_enrich\] to
#' determine the number of APOBEC signature mutations in excess of what would be
#' expected by random mutagenesis. Calculated values are rounded to the nearest
#' whole number.
#' 
#'
#' @param mut_data An data.frame with 5 columns in the following order (column
#' names doesn't matter).
#'  - sample: sample names
#'  - chromosome: The variant affected chromosome (chr1).
#'  - pos: mutation position in the chromosome
#'  - ref: The plus strand reference allele at this position.
#'  - alt: Tumor variant allele at this position.
#' @param ref_genome BSgenome object or name of the installed BSgenome package.
#'  Example: "hg38". Details see "genome" argument of
#'  [getBSgenome][BSgenome::getBSgenome].
#' @param signature_motif A character of the trinucleotide which will be used to
#' calculate MUTATIONS(tCw) and CONTEXT(tCw) to test for enrichment. Defaul is
#' APOBEC signature motif: c("TCA", TCT"). The second base will define
#' background MUTATIONS(CorG) and CONTEXT(CorG). All reverseComplemented motif
#' will also be included.
#' @param background_snv A character specifying the mutation type to restrict
#' the calculation of MUTATIONS(CorG). If NULL, all mutations found in the
#' second base of `signature_motif` will be included. All Complemented
#' substitution will also be included. Default: c("C>T", "C>G").
#' @param bg_extension The number of bases around the variable base to define
#' the background. Default 20L.
#' @param remove_non_background_snv A scalar logical indicates whether mutations
#' not in `background_snv` should be removed before analysis.
#' @return A [data.table][data.table::data.table] describing enrichment per
#' sample.
#' @seealso
#' <https://github.com/PoisonAlien/maftools/blob/master/R/TrinucleotideMatrix.R>
#' @references 
#' - Roberts SA, Lawrence MS, Klimczak LJ, et al. An APOBEC Cytidine Deaminase
#'   Mutagenesis Pattern is Widespread in Human Cancers. Nature genetics.
#'   2013;45(9):970-976. doi:10.1038/ng.2702.
#' - Wang, S., Jia, M., He, Z. et al. APOBEC3B and APOBEC mutational signature
#'   as potential predictive markers for immunotherapy response in non-small
#'   cell lung cancer. Oncogene 37, 3924–3936 (2018).
#'   https://doi.org/10.1038/s41388-018-0245-9 
#' @export
run_motif_fisher <- function(
    mut_data, ref_genome, signature_motif = c("TCA", "TCT"),
    background_snv = c("C>T", "C>G"), remove_non_background_snv = TRUE,
    bg_extension = 20L) {
    assert_pkg("BSgenome")
    assert_pkg("GenomicRanges")
    assert_pkg("GenomeInfoDb")
    assert_pkg("Biostrings")

    assert_class(
        mut_data, function(x) {
            inherits(x, "data.frame") && ncol(x) >= 5L
        },
        msg = "{.cls data.frame} and with at lease 5 columns",
        cross_msg = NULL
    )
    assert_class(
        signature_motif, function(x) {
            is.character(x) && all(nchar(x) == 3L)
        },
        msg = "{.cls character} and all elements must have size {.val 3L}",
        cross_msg = NULL
    )
    assert_class(
        background_snv, function(x) {
            is.character(x) && all(x %chin% names(substitution_pairs))
        },
        msg = "{.cls character} (among {.val {unique(names(substitution_pairs))}})",
        cross_msg = NULL, null_ok = TRUE
    )

    mut_data <- data.table::as.data.table(mut_data)[, .SD, .SDcols = 1:5]
    signature_motif <- union(
        signature_motif,
        as.character(Biostrings::reverseComplement(
            Biostrings::DNAStringSet(signature_motif)
        ))
    )
    signature_mut_bases <- unique(substring(signature_motif, 2L, 2L))
    background_snv_type <- background_snv %||%
        grep(
            sprintf("^(%s)", paste0(signature_mut_bases, collapse = "|")),
            names(substitution_pairs),
            value = TRUE, perl = TRUE
        )
    background_snv_type <- unique(standardize_substitution(background_snv_type))

    # check signature_motif, background_snv, signature_mut_bases are compipable
    signature_motif_mut_bases <- unique(substring(signature_motif, 2L, 2L))
    background_snv_mut_bases <- unique(substring(background_snv_type, 1L, 1L))
    background_snv_mut_bases <- union(
        background_snv_mut_bases,
        as.character(Biostrings::complement(
            Biostrings::DNAStringSet(background_snv_mut_bases)
        ))
    )
    if (!setequal(signature_motif_mut_bases, background_snv_mut_bases)) {
        cli::cli_abort(c(
            "None-compatible bases found",
            i = "{.arg signature_motif} variant bases: {.val {signature_motif_mut_bases}}",
            i = "{.arg background_snv} variant bases: {.val {background_snv_mut_bases}}"
        ))
    }

    ref_genome <- BSgenome::getBSgenome(ref_genome)
    data.table::setnames(
        mut_data, c("sample", "chr", "start", "ref", "alt")
    )
    mut_data[, c("ref", "alt") := lapply(.SD, as.character),
        .SDcols = c("ref", "alt")
    ]
    mut_data <- mut_data[
        !is.na(start) & ref != alt &
            nchar(ref) == 1L & nchar(ref) == 1L
    ]
    if (nrow(mut_data) == 0L) {
        cli::cli_abort("No SNPs to analyze!")
    }
    mut_data[, start := as.numeric(start)]
    mut_data[, end := start]

    data.table::setDF(mut_data)
    mut_gr <- GenomicRanges::makeGRangesFromDataFrame(
        mut_data,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end",
        starts.in.df.are.0based = FALSE
    )
    GenomicRanges::strand(mut_gr) <- "+"

    # match seqlevelstyle between mut_data and ref_genome
    ref_genome_style <- GenomeInfoDb::seqlevelsStyle(ref_genome)
    if (all(GenomeInfoDb::seqlevelsStyle(mut_gr) != ref_genome_style)) {
        cli::cli_inform(
            "Mapping seqnames of {.arg mut_data} to {.arg ref_genome} ({.val {ref_genome_style}})"
        )
        GenomeInfoDb::seqlevelsStyle(mut_gr) <- ref_genome_style
    }
    GenomeInfoDb::seqlevels(mut_gr) <- GenomeInfoDb::seqlevels(ref_genome)
    GenomeInfoDb::seqinfo(mut_gr) <- GenomeInfoDb::seqinfo(ref_genome)

    # Meaure nucleotide frequency and motifs within upstream and downstream of
    # mutated base;
    mut_gr <- mut_gr[GenomicRanges::trim(mut_gr) == mut_gr]

    motif <- suppressWarnings(gr_extend(mut_gr, extension = 1L))
    bg <- suppressWarnings(gr_extend(mut_gr, extension = bg_extension))

    # omit ranges out-of-ranges
    good_ranges <- GenomicRanges::trim(motif) == motif &
        GenomicRanges::trim(bg) == bg
    mut_gr <- mut_gr[good_ranges]
    motif <- BSgenome::getSeq(
        x = ref_genome, motif[good_ranges],
        as.character = FALSE
    )
    bg <- BSgenome::getSeq(
        x = ref_genome, bg[good_ranges],
        as.character = FALSE
    )

    # define substitution and motif (add upstream and downstream bases)
    ## This is nucleotide frequcny and motif frequency
    cli::cli_inform("Defining signature motif frequency and background context")
    compile_data <- cbind(
        sample = mut_gr$sample,
        Biostrings::alphabetFrequency(bg)[, Biostrings::DNA_BASES],
        Biostrings::trinucleotideFrequency(bg),
        define_snv_motif(motif, mut_gr$ref, mut_gr$alt)
    )
    cli::cli_inform(c(
        i = "Using {.val {background_snv_type}} to define {.field background_mut_freq}",
        i = "Using {.val {signature_motif}} to define {.field signature_mut_freq}"
    ))
    compile_data[
        , c("background_mut_freq", "signature_mut_freq") := {
            tmp <- as.character(SubstitutionType) %chin% background_snv_type # nolint
            list(tmp, tmp & as.character(motif) %chin% signature_motif)
        }
    ]
    cli::cli_inform(c(
        i = "Using {.val {signature_mut_bases}} to define {.field background_context}",
        i = "Using {.val {signature_motif}} to define {.field signature_context}"
    ))
    compile_data[, background_context := rowSums(.SD), # nolint
        .SDcols = signature_mut_bases
    ]
    compile_data[, signature_context := rowSums(.SD), # nolint
        .SDcols = signature_motif
    ]
    # only include SNV in the background mutation?
    if (remove_non_background_snv) {
        compile_data <- compile_data[(background_mut_freq)]
    }
    compile_data <- compile_data[,
        c(lapply(.SD, sum), list(n_mutations = .N)),
        .SDcols = !c(
            "Substitution", "SubstitutionMotif", "SubstitutionType",
            "SubstitutionTypeMotif"
        ),
        by = "sample"
    ]
    compile_data[
        , non_signature_mut_freq := n_mutations - signature_mut_freq # nolint
    ]
    compile_data[
        , enrichment := (signature_mut_freq / background_mut_freq) / # nolint
            (signature_context / background_context) # nolint
    ]
    ### One way Fisher test to estimate over representation og APOBEC associated tcw mutations
    cli::cli_inform("Performing one-way Fisher's test")
    compile_data[
        ,
        c("fisher.pvalue", "or", "ci.low", "ci.up") := {
            tmp <- stats::fisher.test(
                matrix(
                    c(
                        .SD$signature_mut_freq,
                        .SD$signature_context,
                        .SD$background_mut_freq - .SD$signature_mut_freq,
                        .SD$background_context - .SD$signature_context
                    ),
                    nrow = 2L
                ),
                alternative = "greater"
            )
            list(tmp$p.value, tmp$estimate, tmp$conf.int[1L], tmp$conf.int[2L])
        },
        by = .I
    ]
    compile_data[, .SD,
        .SDcols = c(
            "sample", "signature_mut_freq", "background_mut_freq",
            "signature_context", "background_context", "non_signature_mut_freq",
            "n_mutations", "enrichment", "fisher.pvalue",
            "or", "ci.low", "ci.up"
        )
    ]
}

gr_extend <- function(x, extension = 1L, use.names = TRUE) {
    GenomicRanges::update_ranges(x,
        start = GenomicRanges::start(x) - extension,
        end = GenomicRanges::end(x) + extension,
        use.names = use.names
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

utils::globalVariables(c(
    "background_context", "signature_context",
    "non_signature_mut_freq", "n_mutations",
    "signature_mut_freq", "background_mut_freq",
    "start", "end", "enrichment", "SubstitutionType",
    "ref", "alt"
))

standardize_substitution <- function(x) {
    unname(substitution_pairs[x])
}

combine_snv_motif <- function(substitution, motif, motif_extension) {
    paste0(
        Biostrings::substring(motif, 1L, motif_extension),
        "[", substitution, "]",
        Biostrings::substring(
            motif, motif_extension + 2L,
            motif_extension * 2L + 1L
        )
    )
}

# values should be the standardized substitution
# names are substitution can be found out
#' @keywords internal
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
