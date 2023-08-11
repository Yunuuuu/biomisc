#' Overdispersion significance analysis for geneset
#' 
#' Principal component analysis was performed on the expression of the specified
#' geneset (`features`), and the variance explained by principal component (PC)1
#' was calculated. For each module, 1,000 random lists of genes with similar
#' expression levels were generated, and the variance explained by PC1 in those
#' genesets was calculated. The significance (p) was calculated as the fraction
#' of random genesets that resulted in a higher PC1 variance than the geneset
#' itself. 
#' @inheritParams sample_sim
#' @inheritDotParams BiocSingular::runPCA -x 
#' @param log A scalar logical indicates whether do logarithms.
#' @references
#' Barkley, D., Moncada, R., Pour, M. et al. Cancer cell states recur across
#' tumor types and form specific interactions with the tumor microenvironment.
#' Nat Genet 54, 1192–1201 (2022). https://doi.org/10.1038/s41588-022-01141-9
#' @seealso
#' <https://github.com/yanailab/PanCancer>
overdispersion <- function(matrix, features, ..., log = FALSE, size = 10^3L, nbin = 25L) {
    assert_pkg("BiocSingular")
    features <- intersect(features, rownames(matrix))
    # https://github.com/yanailab/PanCancer/blob/49e7b270ec55dbe72076a5cae516ff0931fe7fe4/experiments.R#L326
    value <- BiocSingular::runPCA(x = t(matrix[features, , drop = FALSE]),
        rank = 2L, get.pcs = FALSE, get.rotation = FALSE, ...
    )$sdev
    sample_list <- sample_sim(matrix, features, size = size, nbin = nbin)
    samples <- sapply(sample_list, function(sample_features) {
        BiocSingular::runPCA(x = t(matrix[sample_features, , drop = FALSE]),
            rank = 2L, get.pcs = FALSE, get.rotation = FALSE, ...
        )$sdev
    })
    dispersion <- mean(samples[1L, ] >= value[1L])
    coordination <- mean(
        (samples[1L, ] / samples[2L, ]) >= (value[1L] / value[2L])
    )
    out <- c("dispersion" = dispersion, "coordination" = coordination)
    if (isTRUE(log)) {
        out <- log10(out)
    }
    out
}
