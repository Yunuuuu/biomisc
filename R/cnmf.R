#' Consensus Non-negative Matrix factorization
#'
#' A modified version of the cNMF algorithm, implemented in R
#'
#' @param matrix A matrix with row is features and column is samples.
#' @param rank Specification of the factorization rank.
#' @param n_iters The number of times to run NMF internally before making the
#' consensus
#' @param rho Determines number of neighbors to use for calculating KNN distance
#' as rho * n_iters. can be thought of the fraction of replicates that must
#' yield a component approximately matching a program in order for that program
#' to be kept by cNMF.
#' @param min_dist distance threshold that determines how close a component must
#' be to its nearest neighbors in Euclidean space to be considered
#' ‘approximately matching’.
#' @param min_fraction Only featuers with value (above 0) in at least
#' `min_fraction` of samples are used in cNMF.
#' @param silhouette If TRUE, will calculate silhouette score.
#' @inheritDotParams RcppML::nmf -data -k
#' @param cores Parallelization is applied with OpenMP using the number of
#' threads; `0` corresponds to all threads.
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
#' - <https://github.com/dylkot/cNMF>
#' @export
cnmf <- function(matrix, rank, n_iters = 100L, rho = 0.3, min_dist = 0.03, min_fraction = 0.002, silhouette = TRUE, ..., cores = 0L) {
    if (isTRUE(silhouette)) assert_pkg("cluster")
    # RcppML package must be loaded to run nmf
    if (require("RcppML", quietly = TRUE, character.only = TRUE) &&
        utils::packageVersion("RcppML") >= "0.5.5") { # nolint
        cli::cli_abort(sprintf(
            "%s (>=0.5.5) must be installed to use %s.",
            format_pkg("RcppML"), format_fn("cnmf")
        ))
    }
    assert_(rho, function(x) {
        is_scalar_numeric(x) && x > 0L && x <= 1L
    }, "numeric in (0, 1]")

    orig_matrix <- matrix
    matrix <- orig_matrix[
        rowMeans(orig_matrix > 0L) > min_fraction, ,
        drop = FALSE
    ]

    # https://github.com/seanken/CompareSequence/blob/main/ComparePackage_R/CompareSeqR/R/cNMF.R#L53
    cli::cli_inform("Runing NMF")
    old_threads <- options(RcppML.threads = cores)
    on.exit(do.call(`options`, old_threads))
    rank <- as.integer(rank)
    w_list <- lapply(seq_len(n_iters), function(i) {
        RcppML::nmf(data = matrix, k = rank, ...)@w
    })

    cli::cli_inform("Idenfitying consensus programs")
    w <- do.call(cbind, w_list)

    # https://github.com/dylkot/cNMF/blob/master/src/cnmf/cnmf.py
    # Defining consensus w and H
    transposed_w <- t(w) / sqrt(colSums(w^2L))
    L <- as.integer(rho * n_iters)
    # dist regard row as observations
    dist <- stats::dist(transposed_w, method = "euclidean")
    ave_dist <- apply(as.matrix(dist), 1L, function(x) {
        # find the mean over those n_neighbors
        # (excluding self, which has a distance of 0)
        sum(sort(x)[seq_len(L + 1L)]) / L
    }, simplify = TRUE)
    transposed_w <- transposed_w[ave_dist < min_dist, , drop = FALSE]

    # kmeans regard row as observations
    km <- stats::kmeans(transposed_w, centers = rank)
    if (isTRUE(silhouette)) {
        silhouette_score <- cluster::silhouette(
            km$cluster,
            stats::dist(transposed_w, method = "euclidean")
        )
        silhouette_score <- mean(silhouette_score[, 3L, drop = TRUE])
    } else {
        silhouette_score <- NULL
    }
    w_consensus <- data.table::as.data.table(transposed_w)
    w_consensus[, .__groups := km$cluster] # nolint
    w_consensus <- w_consensus[, lapply(.SD, stats::median), by = ".__groups"]
    w_consensus <- as.matrix(w_consensus[, !".__groups"])
    w_consensus <- t(w_consensus / rowSums(abs(w_consensus)))

    h <- RcppML::project(w = w_consensus, data = matrix)
    w_final <- t(RcppML::project(w = h, data = t(orig_matrix)))

    list(w = w_final, h = h, rank = rank, silhouette_score = silhouette_score)
}

utils::globalVariables(".__groups")
