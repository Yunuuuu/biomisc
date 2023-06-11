#' Calculates trinucleotide (or more) context fraction
#'
#' Determines trinucleotide (or more) context fraction from a mutation counts
#' list. 
#' @param mut_data A data.frame of mutation data.
#' @param sample_field A string specifying the sample id column in mut_data.
#' @param chr_field,mut_pos A string, specifying the chromosome and variants
#'  positions column in the mut_data.
#' @param ref_field,alt_field A string specifying the column of reference allele
#'  or variant allele in mut_data.
#' @param extension The number of bases, that's extracted upstream and
#' downstream of the base substitutions. (Default: 1).
#' @param method A string, Default NULL, means no adjustment will be applied.
#' Details see [context_sizefactor]
#' @inheritParams context_sizefactor
#' @return the mutation context matrix
#' @export 
snv_sub_matrix <- function(mut_data, sample_field = NULL, ref_genome = NULL, chr_field = "chr", mut_pos = "pos", ref_field = "ref", alt_field = "alt", extension = 1L, method = NULL, contigs = NULL, exome_ranges = NULL) {
    assert_df_with_columns(mut_data, c(
        sample_field, chr_field, mut_pos, ref_field, alt_field
    ))
    ref_genome <- get_genome(ref_genome)

    # create mutation GenomicRanges
    snv_context_data <- snv_sub_context(
        mut_data,
        chr_field = chr_field,
        mut_pos = mut_pos, ref_field = ref_field,
        alt_field = alt_field,
        ref_genome = ref_genome,
        extension = extension,
        bg_extension = NULL
    )
    snv_sub_context <- combine_sub_context(
        snv_context_data$sub_type_context,
        sub = snv_context_data$sub_type,
        extension = extension
    )
    standard_snv_sub_context <- enumerate_standard_snv_sub_context(extension)
    snv_sub_type_context <- factor(snv_sub_context, standard_snv_sub_context)

    if (is.null(sample_field)) {
        sub_type_matrix <- c(table(snv_sub_type_context))
        sub_type_matrix <- matrix(sub_type_matrix, ncol = 1L)
        rownames(sub_type_matrix) <- standard_snv_sub_context
    } else {
        sub_type_matrix <- table(
            snv_sub_type_context,
            S4Vectors::mcols(snv_context_data$mut)[[sample_field]]
        )
    }

    if (!is.null(method)) {
        snv_context <- sub_context_to_mut_context(standard_snv_sub_context)
        sub_type_matrix <- sub_type_matrix *
            calculate_context_sizefactor(snv_context,
                method = method, contigs = contigs,
                ref_genome = ref_genome, exome_ranges = exome_ranges
            )
    }
    sub_type_matrix
}
