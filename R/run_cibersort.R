#' a modified CIBERSORT function
#' @description  Robust enumeration of cell subsets from tissue expression
#'   profiles
#' @details By default, CIBERSORT estimates the relative fraction of each cell
#'   type in the signature matrix, such that the sum of all fractions is equal
#'   to 1 for a given mixture sample.
#'
#'   CIBERSORT can also be used to produce a score that quantitatively measures
#'   the overall abundance of each cell type (as described in
#'   \href{https://www.nature.com/articles/nmeth.3337}{Analysis of deconvolution
#'   consistency}). Briefly, the absolute immune fraction score is estimated by
#'   the median expression level of all genes in the signature matrix divided by
#'   the median expression level of all genes in the mixture. more details see
#'   \href{https://cibersort.stanford.edu/manual.php#run}{CIBERSORT Manual}
#'
#' @param mixture_data a `matrix` object whose rownames is the gene names and is
#' consistent with `sig_data` gene names; Data should be in non-log space.
#' Note: if maximum expression value is <50; CIBERSORT will assume that data are
#' in log space, and will anti-log all expression values by 2x. If gene symbols
#' are redundant, CIBERSORT will choose the one with highest mean expression
#' across the mixtures. CIBERSORT performs a feature selection and therefore
#' typically does not use all genes in the signature matrix. It is generally ok
#' if some genes are missing from the user’s mixture file. If <50% of signature
#' matrix genes overlap, CIBERSORT will issue a warning.  Normal quantification
#' of RNA-seq like FPKM, and TPM can be used.
#' @param sig_data CIBERSORT requires an input matrix of reference gene
#' expression signatures, or signature matrix with gene names in rownames, for
#' routine analysis. This is stored in a Signature Genes File and consists of a
#' table with groups of "barcode" genes whose expression values collectively
#' define a unique gene expression signature for each component pure cell
#' population that will be used to deconvolute the mixture. if \code{NULL},
#' \code{LM22} with HUGO gene symbols as gene names will be used.
#' @param perm Number of permutations; set to >=100 to calculate p-values.
#'   default: `200`.
#' @param quantile_norm Quantile normalization of input mixture. default: `TRUE`
#' @param absolute Run CIBERSORT in absolute mode default: \code{FALSE}. note
#'   that cell subsets will be scaled by their absolute levels and will not be
#'   represented as fractions (to derive the default output, normalize absolute
#'   levels such that they sum to 1 for each mixture sample); the sum of all
#'   cell subsets in each mixture sample will be added to the ouput ('Absolute
#'   score'). If LM22 is used, this score will capture total immune content.
#' @param abs_method if absolute is set to \code{TRUE}, abs_method choose
#'   method: 'no_sumto1' or 'sig_score' - sig_score = for each mixture sample,
#'   define S as the median expression level of all genes in the signature
#'   matrix divided by the median expression level of all genes in the mixture.
#'   Multiple cell subset fractions by S. - no_sumto1 = remove sum to 1
#'   constraint
#' @return a data.frame
#' @references
#'  - Newman, A., Liu, C., Green, M. et al. Robust enumeration of cell subsets
#'   from tissue expression profiles. Nat Methods 12, 453–457 (2015).
#'   <https://doi.org/10.1038/nmeth.3337>
#'  - Chen B., Khodadoust M.S., Liu C.L., Newman A.M., Alizadeh A.A. (2018)
#'   Profiling Tumor Infiltrating Immune Cells with CIBERSORT. In: von Stechow
#'   L. (eds) Cancer Systems Biology. Methods in Molecular Biology, vol 1711.
#'   Humana Press, New York, NY.  <https://doi.org/10.1007/978-1-4939-7493-1_12>
#' @export
run_cibersort <- function(mixture_data, sig_data = NULL,
                          perm = 200L, quantile_norm = TRUE, absolute = FALSE,
                          abs_method = "sig_score") {
    assert_pkg("e1071")
    assert_pkg("preprocessCore")

    if (absolute) {
        abs_method <- match.arg(abs_method, c("no_sumto1", "sig_score"))
    }

    # choose the default signature data or provided by users

    if (is.null(sig_data)) {
        sig_data <- cibersort_lm22
    }

    if (max(mixture_data) < 50L) {
        cli::cli_warn(
            "find {.arg mixture_data} is logged, we'll anti-log by {.code 2^mixture_data}!"
        )
        mixture_data <- 2^mixture_data
    }

    # quantile normalization of mixture file
    if (quantile_norm) {
        tmpc <- colnames(mixture_data)
        tmpr <- rownames(mixture_data)
        mixture_data <- preprocessCore::normalize.quantiles(mixture_data)
        colnames(mixture_data) <- tmpc
        rownames(mixture_data) <- tmpr
    }

    # store original mixtures median value
    mixture_data_median <- max(stats::median(mixture_data), 1L)

    # intersect genes and keep them in the same order
    common_genes <- intersect(row.names(sig_data), row.names(mixture_data))
    sig_data <- sig_data[common_genes, ]
    mixture_data <- mixture_data[common_genes, ]

    # standardize sig matrix
    sig_data <- (sig_data - mean(sig_data)) / stats::sd(sig_data)

    # empirical null distribution of correlation coefficients
    if (perm > 0L) {
        cli::cli_inform("Calculating null distribution by permutation...")
        nulldist_ecdf <- stats::ecdf(
            cibersort_do_perm(
                perm, sig_data, mixture_data,
                absolute, abs_method
            )
        )
    }
    cli::cli_inform("Running CIBERSORT analysis...")
    p <- progressr::progressor(steps = ncol(mixture_data))
    results <- future.apply::future_apply(mixture_data, 2L, function(sample_col) {
        p()
        # standardize mixture
        scale_sample_col <- scale(sample_col, center = TRUE, scale = TRUE)

        # run SVR core algorithm
        result <- cibersort_core_algorithm(
            sig_data = sig_data, y = scale_sample_col,
            absolute, abs_method
        )

        if (absolute && identical(abs_method, "sig_score")) {
            result$w <- result$w * stats::median(sample_col) /
                mixture_data_median
        }

        # calculate p-value
        if (perm > 0L) {
            pval <- 1L - nulldist_ecdf(result$mix_r)
        } else {
            pval <- 1L
        }

        # return results
        if (absolute) {
            c(result$w, pval, result$mix_r, result$mix_rmse, sum(result$w))
        } else {
            c(result$w, pval, result$mix_r, result$mix_rmse)
        }
    }, simplify = FALSE, future.globals = TRUE)

    # return matrix object containing all results
    results <- do.call("rbind", results)
    rownames(results) <- colnames(mixture_data)
    if (!absolute) {
        colnames(results) <- c(
            colnames(sig_data),
            "P-value", "Correlation", "RMSE"
        )
    } else {
        colnames(results) <- c(
            colnames(sig_data),
            "P-value", "Correlation", "RMSE",
            paste("Absolute score (", abs_method, ")", sep = "")
        )
    }
    results <- data.table::as.data.table(results, keep.rownames = "Samples")
    data.table::setDF(results)
    results
}

# Core algorithm
cibersort_core_algorithm <- function(sig_data, y, absolute, abs_method) {
    out <- lapply(c(0.25, 0.5, 0.75), function(nu) {
        e1071::svm(
            sig_data, y,
            type = "nu-regression",
            kernel = "linear", nu = nu,
            scale = FALSE
        )
    })
    corrv <- nusvm <- numeric(length(out))

    # do cibersort
    for (i in seq_along(out)) {
        weights <- t(out[[i]]$coefs) %*% out[[i]]$SV
        weights[which(weights < 0L)] <- 0L
        w <- weights / sum(weights)
        u <- sweep(sig_data, MARGIN = 2L, w, "*")
        k <- apply(u, 1L, sum)
        nusvm[[i]] <- sqrt((mean((k - y)^2L)))
        corrv[[i]] <- stats::cor(k, y)
    }

    # pick the best model (nusvm is the rmse)
    mn <- which.min(nusvm)
    model <- out[[mn]]

    # get and normalize coefficients
    q <- t(model$coefs) %*% model$SV
    q[which(q < 0L)] <- 0L
    # relative space (returns fractions)
    if (!absolute || identical(abs_method, "sig_score")) w <- (q / sum(q))
    # absolute space (returns scores)
    if (absolute && identical(abs_method, "no_sumto1")) w <- q
    list(w = w, mix_rmse = nusvm[[mn]], mix_r = corrv[[mn]])
}

# do permutations
cibersort_do_perm <- function(perm, sig_data, mixture_data, absolute, abs_method) {
    p <- progressr::progressor(steps = perm)
    future.apply::future_vapply(
        seq_len(perm), function(i) {
            p(message = sprintf("Permuatating %d times", i))
            # random mixture
            yr <- mixture_data[sample.int(length(mixture_data), nrow(sig_data))]

            # standardize mixture
            yr <- scale(yr, center = TRUE, scale = TRUE)

            # run CIBERSORT core algorithm
            cibersort_core_algorithm(sig_data, yr, absolute, abs_method)$mix_r
        }, numeric(1L),
        USE.NAMES = FALSE,
        future.seed = TRUE, future.globals = TRUE
    )
}
