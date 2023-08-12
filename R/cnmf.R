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
#' @param rho Determines number of neighbors to use for calculating KNN distance
#' as rho * n_iters. can be thought of the fraction of replicates that must
#' yield a component approximately matching a program in order for that program
#' to be kept by cNMF.
#' @param min_dist distance threshold that determines how close a component must
#' be to its nearest neighbors in Euclidean space to be considered
#' ‘approximately matching’.
#' @return A list.
#' @references
#' - Simmons, S.K., Lithwick-Yanai, G., Adiconis, X. et al. Mostly natural
#'   sequencing-by-synthesis for scRNA-seq using Ultima sequencing. Nat
#'   Biotechnol 41, 204–211 (2023). https://doi.org/10.1038/s41587-022-01452-6
#' - Dylan KotliarAdrian VeresM Aurel NagyShervin TabriziEran HodisDouglas A
#'   MeltonPardis C Sabeti (2019) Identifying gene expression programs of
#'   cell-type identity and cellular activity with single-cell RNA-Seq eLife
#'   8:e43803. https://doi.org/10.7554/eLife.43803
#' @seealso
#' - <https://github.com/seanken/CompareSequence>
#' - <https://github.com/dylkot/cNMF/blob/master/src/cnmf/cnmf.py>
#' @export
cnmf <- function(matrix, min_fraction = 0.002, k = 15L, n_iters = 100L, rho = 0.3, min_dist = 0.03) {
    assert_pkg("RcppML")
    assert_pkg("cluster")
    assert_class(matrix, is.matrix, "{.cls matrix}")
    assert_class(rho, function(x) {
        is.numeric(x) && x > 0L && x <= 1L
    }, "(0, 1] {.cls numeric}", cross_msg = NULL)

    orig_matrix <- matrix
    matrix <- orig_matrix[
        rowMeans(orig_matrix > 0L) > min_fraction, ,
        drop = FALSE
    ]
    # https://github.com/seanken/CompareSequence/blob/main/ComparePackage_R/CompareSeqR/R/cNMF.R#L53
    cli::cli_inform("Runing NMF")
    w_list <- lapply(seq_len(n_iters), function(i) {
        RcppML::nmf(A = matrix, k = k)$w
    })
    
    cli::cli_inform("Idenfity consensus programs")
    w <- do.call(cbind, w_list)

    # Defining consensus w and H
    transposed_w <- t(w) / sqrt(colSums(w^2L))
    L <- as.integer(rho * n_iters)
    # dist regard row as observations
    dist <- stats::dist(transposed_w, method = "euclidean")
    ave_dist <- apply(as.matrix(dist), 1L, function(x) {
        # find the mean over those n_neighbors
        # (excluding self, which has a distance of 0)
        sum(sort(x, , decreasing = FALSE)[seq_len(L + 1L)]) / L
    }, simplify = TRUE)
    transposed_w <- transposed_w[ave_dist < min_dist, , drop = FALSE]

    # kmeans regard row as observations
    km <- stats::kmeans(transposed_w, centers = k)
    silhouette_score <- cluster::silhouette(km$cluster, dist)
    silhouette_score <- mean(silhouette_score[, 3L, drop = TRUE])
    w_consensus <- data.table::as.data.table(transposed_w)
    w_consensus[, .__groups := km$cluster] # nolint
    w_consensus <- w_consensus[, lapply(.SD, median), by = ".__groups"]
    w_consensus <- as.matrix(w_consensus[, !".__groups"])
    w_consensus <- t(w_consensus / rowSums(abs(w_consensus)))

    h <- RcppML::project(matrix, w = w_consensus)

    w_final <- RcppML::project(orig_matrix, h = h)

    list(w = t(w_final), h = t(h), silhouette_score = silhouette_score)
}
