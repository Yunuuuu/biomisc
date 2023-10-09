#' Inferring Correlation Networks from Genomic Survey Data
#' @description A reimplementation of SparCC algorithm (Friedman et Alm 2012,
#' PLoS Comp Bio, 2012). SparCC is mainly used for calculating correlations in
#' compositional data.
#' @inheritParams SpiecEasi::sparcc
#' @inheritDotParams SpiecEasi::sparcc -data
#' @param times Number of bootstraps
#' @return A list.
#' @seealso
#' [sparcc][SpiecEasi::sparcc]
#' @references 
#' <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687>
#' @export
sparcc <- function(data, ..., times = 200L) {
    assert_pkg("SpiecEasi")
    sparcc0 <- SpiecEasi::sparcc(data = data, ...)

    # bootstrap to calculate P-value ----------------------
    t0 <- SpiecEasi::triu(sparcc0$Cor)
    n <- NROW(data)
    bootstrap_list <- future.apply::future_lapply(seq_len(times), function(i) {
        out <- SpiecEasi::sparcc(
            data = data[sample(n, size = n, replace = TRUE), , drop = FALSE],
            ...
        )
        SpiecEasi::triu(out$Cor)
    }, future.seed = TRUE)
    bootstrap_out <- do.call("rbind", bootstrap_list)
    permutation_list <- future.apply::future_lapply(seq_len(times), function(i) {
        out <- SpiecEasi::sparcc(data = apply(data, 2L, sample), ...)
        SpiecEasi::triu(out$Cor)
    }, future.seed = TRUE)
    permutation_out <- do.call("rbind", permutation_list)

    tmeans <- colMeans(permutation_out)
    ind95 <- max(1, round(0.025 * times)):round(0.975 * times)
    boot_ord <- apply(bootstrap_out, 2L, sort)
    boot_ord95 <- boot_ord[ind95, , drop = FALSE]
    outofrange <- vapply(seq_along(t0), function(i) {
        aitvar <- t0[i]
        range <- range(boot_ord95[, i, drop = TRUE])
        range[1L] > aitvar || range[2L] < aitvar
    }, logical(1L))
    bs_above <- vapply(seq_len(ncol(bootstrap_out)), function(i) {
        sum(bootstrap_out[, i, drop = TRUE] > tmeans[i], na.rm = TRUE)
    }, integer(1L))
    pvals <- ifelse(bs_above > times / 2L,
        2 * (1 - bs_above / times), 2 * bs_above / times
    )
    pvals[pvals > 1] <- 1
    pvals[outofrange] <- NaN

    # coerce pvals into a matrix ----------------------
    pvals_mat <- diag(rep_len(0L, nrow(data)))
    pvals_mat[upper.tri(pvals_mat, diag = FALSE)] <- pvals
    pvals_mat <- pvals_mat + t(pvals_mat)
    diag(pvals_mat) <- 0L

    # return SparCC and P-value
    colnames(sparcc0$Cor) <- rownames(sparcc0$Cor) <- rownames(data)
    colnames(pvals_mat) <- rownames(pvals_mat) <- rownames(data)
    c(sparcc0, list(pvals = pvals_mat))
}
