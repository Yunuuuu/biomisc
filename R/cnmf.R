#' Modified cNMF
#'
#' A modified version of the cNMF algorithm, implemented in R
#'
#' @param matrix A matrix with row is features and column is samples.
#' @param min_fraction Only featuers with value (above 0) in at least
#' `min_fraction` of samples are used in cNMF.
#' @param k Number of programs.
#' @param n_iters The number of times to run NMF internally before making the
#' consensus
#' @param local_neighborhood_size Determines number of neighbors to use for
#' calculating KNN distance as local_neighborhood_size * n_iters.
#' @param min_dist The distance used for filtering out outliers
#' @param numb_smp Number of cells to select, can be a floating number between 0
#' and 1, which means the fraction of cells to choose (if NULL or <=0 no
#' downsampling will be applied).
#' @return A list.
#' @references
#' Simmons, S.K., Lithwick-Yanai, G., Adiconis, X. et al. Mostly natural
#' sequencing-by-synthesis for scRNA-seq using Ultima sequencing. Nat Biotechnol
#' 41, 204–211 (2023). https://doi.org/10.1038/s41587-022-01452-6
#' @seealso 
#' - <https://github.com/seanken/CompareSequence>
#' - <https://github.com/dylkot/cNMF/blob/master/src/cnmf/cnmf.py>
#' @export
cnmf <- function(matrix, min_fraction = 0.002, k = 15L, n_iters = 100L, local_neighborhood_size = 0.3, min_dist = 0.03, numb_smp = NULL) {
    assert_class(matrix, is.matrix, "{.cls matrix}")
    if (!is.null(numb_smp) && numb_smp > 0L) {
        if (numb_smp <= 1L) {
            numb_smp <- ncol(matrix) * numb_smp
        }
        numb_smp <- as.integer(numb_smp)
        numb_smp <- min(numb_smp, ncol(matrix))
        cli::cli_inform("Downsampling {numb_smp} sample{?s}")
        orig_matrix <- matrix[, sample(ncol(matrix), numb_smp), drop = FALSE]
    } else {
        orig_matrix <- matrix
    }
    matrix <- orig_matrix[
        rowMeans(orig_matrix > 0L) > min_fraction, ,
        drop = FALSE
    ]
    # https://github.com/seanken/CompareSequence/blob/main/ComparePackage_R/CompareSeqR/R/cNMF.R#L53
    cli::cli_inform("Runing NMF")
    w_list <- lapply(seq_len(n_iters), function(i) {
        RcppML::nmf(matrix, k)$w
    })

    cli::cli_inform("Combining")
    W <- do.call(cbind, w_list)
    W <- t(t(W) / sqrt(colSums(W^2L)))

    dist <- 2L - 2L * (t(W) %*% W)
    L <- local_neighborhood_size * n_iters
    ave_dist <- apply(dist, 1L, function(x) {
        mean(sort(x, , decreasing = FALSE)[seq_len(L)])
    })
    W <- t(W[, ave_dist < min_dist])

    factors <- stats::kmeans(W, k)

    W_new <- data.table::as.data.table(W)
    W_new[, .__clusters := factors$cluster] # nolint
    W_new <- W_new[, lapply(.SD, median), by = ".__clusters"]
    W_new <- as.matrix(W_new[, !".__clusters"])
    W_new <- t(W_new / rowSums(abs(W_new)))

    H <- RcppML::project(matrix, w = W_new)

    W_fin <- RcppML::project(orig_matrix, h = H)

    list(H = t(H), W = t(W_fin))
}
