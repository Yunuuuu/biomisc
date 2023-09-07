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

    assert_(
        mut_data, function(x) {
            inherits(x, "data.frame") && ncol(x) >= 5L
        },
        what = "a {.cls data.frame} (with at lease 5 columns)"
    )
    assert_(
        signature_motif, function(x) {
            is.character(x) && all(nchar(x) == 3L)
        },
        what = "an atomic {.cls character} (all have nchar 3L)"
    )
    assert_(
        background_snv, function(x) {
            is.character(x) && all(x %chin% names(standard_snv_sub_pairs))
        },
        what = sprintf(
            "an atomic {.cls character} (all must in %s})",
            oxford_comma(unique(names(standard_snv_sub_pairs)))
        ),
        null_ok = TRUE
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
            names(standard_snv_sub_pairs),
            value = TRUE, perl = TRUE
        )
    background_snv_type <- unique(standardize_snv_sub(background_snv_type))

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
    mut_gr <- map_seqnames(mut_gr, ref_genome_style, "mut_data")
    GenomeInfoDb::seqlevels(mut_gr) <- GenomeInfoDb::seqlevels(ref_genome)
    suppressWarnings(
        GenomeInfoDb::seqinfo(mut_gr) <- GenomeInfoDb::seqinfo(ref_genome)
    )

    # Meaure nucleotide frequency and motifs within upstream and downstream of
    # mutated base;
    out_ranges <- GenomicRanges::trim(mut_gr) != mut_gr
    if (any(out_ranges)) {
        cli::cli_warn("Removing {.val {sum(out_ranges)}} ({format_percent(mean(out_ranges))}) out-of-bound mutations compared to {.arg ref_genome}")
    }
    mut_gr <- mut_gr[!out_ranges]

    # Get motif and background base sequencing
    motif <- suppressWarnings(granges_extend(mut_gr, extension = 1L))
    bg <- suppressWarnings(granges_extend(mut_gr, extension = bg_extension))

    # omit out-of-bound ranges
    good_ranges <- GenomicRanges::trim(bg) == bg

    mut_gr <- mut_gr[good_ranges]
    motif <- BSgenome::getSeq(
        x = ref_genome, motif[good_ranges],
        as.character = FALSE
    )
    # check mutation reference allele match the second base in motif
    matched_ranges <- mut_gr$ref == as.character(
        Biostrings::subseq(motif, 2L, 2L)
    )

    if (any(!matched_ranges)) {
        cli::cli_warn("Removing {.val {sum(!matched_ranges)}} ({format_percent(mean(!matched_ranges))}) mutations whose reference allele cannot match the {.arg ref_genome}")
    }

    mut_gr <- mut_gr[matched_ranges]
    motif <- motif[matched_ranges]
    bg <- BSgenome::getSeq(
        x = ref_genome,
        bg[good_ranges][matched_ranges],
        as.character = FALSE
    )

    # define substitution and motif (add upstream and downstream bases)
    ## This is nucleotide frequcny and motif frequency
    cli::cli_inform("Defining signature motif frequency and background context")
    compile_data <- cbind(
        sample = mut_gr$sample,
        Biostrings::alphabetFrequency(bg)[, Biostrings::DNA_BASES],
        Biostrings::trinucleotideFrequency(bg),
        define_snv_sub_tri_context(motif, mut_gr$ref, mut_gr$alt)
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

utils::globalVariables(c(
    "background_context", "signature_context",
    "non_signature_mut_freq", "n_mutations",
    "signature_mut_freq", "background_mut_freq",
    "start", "end", "enrichment", "SubstitutionType",
    "ref", "alt"
))
