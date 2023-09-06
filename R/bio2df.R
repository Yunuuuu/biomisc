#' Transfrom a bioconductor Object into a data.frame
#' 
#' It's often we want to combine the expression data and phenotypic (or feature)
#' data into a data.frame to do something like plot (ggplot2) or simple
#' statistical test.
#' @param x A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment]
#'   or [ExpressionSet][Biobase::ExpressionSet] object.
#' @param ... For `bio2df`, other arguments passed to specific methods.
#' @name bio2df
#' @return A data.frame.
#' @export
bio2df <- function(x, ...) {
    UseMethod("bio2df")
}

#' @param assay Integer scalar or string indicating which assay of `x` to use.
#' @param row_as_observations If TRUE, a data.frame combined the assay data and
#'   `colData(x)` will be returned, otherwise, a data.frame combined the
#'   transposed assay data and `rowData(x)`.
#' @rdname bio2df
#' @export
bio2df.SummarizedExperiment <- function(x, ..., assay = NULL, row_as_observations = FALSE) {
    assay <- assay %||% 1L
    assay_data <- SummarizedExperiment::assay(x, i = assay)
    if (isTRUE(row_as_observations)) {
        assay_data <- data.frame(assay_data,
            check.names = FALSE, fix.empty.names = FALSE
        )
        feature_data <- SummarizedExperiment::rowData(x)
        feature_data <- data.frame(
            feature_data,
            check.names = FALSE, fix.empty.names = FALSE
        )
        cbind(assay_data, feature_data)
    } else {
        assay_data <- data.frame(t(assay_data),
            check.names = FALSE, fix.empty.names = FALSE
        )
        pheno_data <- SummarizedExperiment::colData(x)
        pheno_data <- data.frame(
            pheno_data,
            check.names = FALSE, fix.empty.names = FALSE
        )
        cbind(assay_data, pheno_data)
    }
}
