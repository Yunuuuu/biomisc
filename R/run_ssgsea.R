#' Single-sample GSEA
#'
#' Project each sample within a data set onto a space of gene set enrichment
#' scores using the ssGSEA projection methodology. According to GenePattern
#' (modified) and GSVA
#'
#' @param data_expr gene expression data. For ssGSEA, normalizd exression values
#'   with gene length adjusted were needed.
#' @param gene_set_list gene sets can be provided as `list` object
#' @param NES whether to calculate normalized enrichment scores
#' @param perm the number of permutations to calculate NES
#' @param sample_norm_type Normalization method applied to expression data.
#'   Supported methods are rank, log.rank, and log.  (Default: rank)
#' @param weight Exponential weight employed in calculation of enrichment
#'   scores. The default value of \code{0.75} was selected after extensive
#'   testing. The module authors strongly recommend against changing from
#'   default. (Default: `0.75`)
#' @param min_sz minimal size of geneSet for testing. (Default: 1L)
#' @param max_sz maximal size of geneSet for testing. (Default: Inf)
#' @return a `list` of projection results for each gene set and each sample
#' @references \itemize{\item Subramanian A, Tamayo P, Mootha VK, Mukherjee S,
#'   Ebert BL, Gillette MA, Paulovich A, Pomeroy SL, Golub TR, Lander ES,
#'   Mesirov JP. Gene set enrichment analysis: A knowledge-based approach for
#'   interpreting genome-wide expression profiles. PNAS.
#'   2005;102(43):15545-15550. \url{https://doi.org/10.1073/pnas.0506580102}
#'   \item Barbie, D., Tamayo, P., Boehm, J. et al. Systematic RNA interference
#'   reveals that oncogenic KRAS-driven cancers require TBK1. Nature 462,
#'   108–112 (2009). \url{https://doi.org/10.1038/nature08460}}
#' @export 
run_ssgsea <- function(data_expr, gene_set_list, NES = TRUE, perm = 1000L,
                       # normalization method applied to input feature data:
                       # "none", "rank", "log" or "log.rank"
                       sample_norm_type = c("rank", "log", "log.rank", "none"),
                       # exponential weight applied to ranking in calculation of
                       # enrichment score
                       weight = 0.75, min_sz = 1L, max_sz = Inf) {
    ## arguments-------------------------------------------

    sample_norm_type <- match.arg(sample_norm_type)

    ## keep gene and sample names -------------------------------------------

    keep_gene_name <- rownames(data_expr)
    keep_sample_name <- colnames(data_expr)

    #############################################
    # Sample normalization
    #############################################

    if (identical(sample_norm_type, "none")) {
        rlang::inform("No normalization would be made")
    } else if (identical(sample_norm_type, "rank")) {
        data_expr <- apply(data_expr, 2L, function(x) {
            as.integer(rank(x, ties.method = "average"))
        })
        data_expr <- 10000L * data_expr / nrow(data_expr)
    } else if (identical(sample_norm_type, "log.rank")) {
        data_expr <- apply(data_expr, 2L, function(x) {
            as.integer(rank(x, ties.method = "average"))
        })
        data_expr <- log(10000L * data_expr / nrow(data_expr) + exp(1L))
    } else if (identical(sample_norm_type, "log")) {
        data_expr[data_expr < 1L] <- 1L
        data_expr <- log(data_expr + exp(1L))
    }

    ## recover gene and sample names -------------------------------------

    rownames(data_expr) <- keep_gene_name
    colnames(data_expr) <- keep_sample_name

    # prune gene sets - only keep genes in data_expr
    gene_set_list <- lapply(gene_set_list, intersect, keep_gene_name)

    ## remove gene sets from the analysis for which no features are available
    ## and meet the minimum and maximum gene-set size specified by the user
    gene_set_list <- filter_genesets(
        gene_set_list,
        min_sz = min_sz,
        max_sz = max_sz
    )

    ## run ssGSEA analysis --------------------------------------------------
    # projecting expression data to gene_set_list
    # data.matrix containing gene expression data
    # exponential weight applied to ranking in calculation of enrichment score
    #################################################

    project_to_geneset(
        data_matrix = data_expr,
        gene_set_list = gene_set_list,
        weight = weight, NES = NES, perm = perm
    )
}

# utility functions -------------------------------------------------------

project_to_geneset <- function(data_matrix, gene_set_list, weight, NES, perm) {
    gene_name <- rownames(data_matrix)
    rlang::inform("Running ssGSEA.............")
    p <- progressr::progressor(steps = ncol(data_matrix))
    es_res <- future.apply::future_lapply(
        seq_len(ncol(data_matrix)), function(j) {
            p(message = "Estimating ssGSEA sccores")
            gene_list <- data_matrix[, j, drop = TRUE]
            names(gene_list) <- gene_name
            gene_list <- sort(gene_list, decreasing = TRUE)

            enrich_score <- lapply(gene_set_list, function(gene_set) {
                if (NES) {
                    ssgsea_nes(
                        gene_list = gene_list,
                        gene_set = gene_set,
                        weight = weight,
                        perm = perm
                    )
                } else {
                    c(ES = ssgsea(
                        gene_list = gene_list,
                        gene_set = gene_set,
                        weight = weight
                    ))
                }
            })
            # when NES is TRUE, a list of three elements: ES, NES, pvalue
            # Otherwise, a list of one element: ES
            data.table::transpose(enrich_score)
        },
        future.globals = TRUE
    )
    ES <- lapply(es_res, "[[", 1L)
    ES <- do.call("cbind", ES)
    rownames(ES) <- names(gene_set_list)
    colnames(ES) <- colnames(data_matrix)
    res <- list(ES = ES)

    if (NES) {
        NES <- lapply(es_res, "[[", 2L)
        NES <- do.call("cbind", NES)
        rownames(NES) <- names(gene_set_list)
        colnames(NES) <- colnames(data_matrix)

        pvalue <- lapply(es_res, "[[", 3L)
        pvalue <- do.call("cbind", pvalue)
        rownames(pvalue) <- names(gene_set_list)
        colnames(pvalue) <- colnames(data_matrix)

        res <- c(res, list(NES = NES, pvalue = pvalue))
    }
    return(res)
}

## optimized version of the function .rndWalk by Alexey Sergushichev
## https://github.com/rcastelo/GSVA/pull/15
## based on his paper https://doi.org/10.1101/060012
ssgsea <- function(gene_list, gene_set, weight) {
    n <- length(gene_list)
    k <- length(gene_set)
    idx <- match(gene_set, names(gene_list))
    gene_list <- abs(gene_list)^weight
    stepCDFinGeneSet <- sum(gene_list[idx] * (n - idx + 1)) /
        sum(gene_list[idx])
    stepCDFoutGeneSet <- (n * (n + 1) / 2 - sum(n - idx + 1)) / (n - k)
    stepCDFinGeneSet - stepCDFoutGeneSet
}

ssgsea_nes <- function(gene_list, gene_set, weight, perm) {
    es <- ssgsea(
        gene_list = gene_list,
        gene_set = gene_set,
        weight = weight
    )
    perm_scores <- vapply(seq_len(perm), function(i) {
        perm_ssgsea_es(
            gene_list = gene_list,
            gene_set = gene_set,
            weight = weight
        )
    }, numeric(1L))

    if (es >= 0L) {
        m <- mean(perm_scores[perm_scores >= 0L])
    } else {
        m <- abs(mean(perm_scores[perm_scores < 0L]))
    }
    NES <- es / m

    if (is.na(NES)) {
        pvalue <- NA_real_
    } else if (NES >= 0L) {
        pvalue <- (sum(perm_scores >= NES) + 1L) / (sum(perm_scores >= 0L) + 1L)
    } else {
        pvalue <- (sum(perm_scores <= NES) + 1L) / (sum(perm_scores < 0L) + 1L)
    }
    return(c(ES = es, NES = NES, pvalue = pvalue))
}

perm_ssgsea_es <- function(gene_list, gene_set, weight) {
    names(gene_list) <- names(gene_list)[
        sample.int(length(gene_list), replace = FALSE)
    ]
    ssgsea(
        gene_list = gene_list,
        gene_set = gene_set,
        weight = weight
    )
}
filter_genesets <- function(genesets, min_sz, max_sz) {
    genesets_sz <- lengths(genesets)
    genesets[
        genesets_sz >= min_sz & genesets_sz <= max_sz
    ]
}
