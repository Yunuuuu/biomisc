#' Sample features with similar expression level
#' @param matrix Normalized gene expression matrix with row is genes and column
#' is samples. 
#' @param features The template for gene sampling, considering genes with
#' comparable expression levels of these genes. 
#' @param size A non-negative integer giving the number of sampling.
#' @param nbin Divides the dataset into approximately `nbin` clusters based on
#' gene expression levels. Subsequently, the sampling process will be applied
#' within each cluster for the specified features.
#' @return A list of character for each sampled genes.
#' @references
#' Barkley, D., Moncada, R., Pour, M. et al. Cancer cell states recur across
#' tumor types and form specific interactions with the tumor microenvironment.
#' Nat Genet 54, 1192–1201 (2022). https://doi.org/10.1038/s41588-022-01141-9
#' @seealso
#' <https://github.com/yanailab/PanCancer>
#' @export
sample_sim <- function(matrix, features, size = 10^3L, nbin = 25L) {
    assert_class(matrix, is.matrix, "{.cls matrix}")
    features <- intersect(features, rownames(matrix))
    data_avg <- rowMeans(matrix)
    if (nbin >= nrow(matrix)) {
        nbin <- as.integer(nrow(matrix) / 3L)
    }

    # <https://github.com/yanailab/PanCancer/blob/49e7b270ec55dbe72076a5cae516ff0931fe7fe4/seurat_functions_public.R#L1500> 
    data_cut <- cut_number(x = data_avg, n = nbin, labels = FALSE)
    names(data_cut) <- names(data_avg)
    data_group <- split(names(data_cut), data_cut)
    feature_group <- split(features, data_cut[features])
    idx_group <- split(seq_along(features), data_cut[features])
    lapply(seq_len(size), function(i) {
        for (id in names(feature_group)) {
            features[idx_group[[id]]] <- sample(
                data_group[[id]],
                size = length(feature_group[[id]])
            )
        }
        features
    })
}

# copied from ggplot2
cut_number <- function(x, n = NULL, ...) {
    probs <- seq(0, 1, length.out = n + 1)
    brk <- stats::quantile(x, probs, na.rm = TRUE)
    if (anyDuplicated(brk)) {
        cli::cli_abort("Insufficient data values to produce {n} bins.")
    }
    cut(x, brk, include.lowest = TRUE, ...)
}
