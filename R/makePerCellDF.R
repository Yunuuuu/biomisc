#' Create a per-cell data.frame
#'
#' Create a per-cell data.frame (i.e., where each row represents a sample /
#' cell) from a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment], most
#' typically for creating custom ggplot2 plots. 
#' @param x A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment]
#'   or [ExpressionSet][Biobase::ExpressionSet] object.
#' @param ... Other arguments passed to specific methods.
#' @return A data.frame.
#' @name makePerCellDF
#' @export
methods::setGeneric(
    "makePerCellDF",
    function(x, ...) standardGeneric("makePerCellDF"),
    signature = "x"
)

#' @param features Character or integer vector specifying the features for which
#' to extract expression profiles across samples. 
#' @param assay.type String or integer scalar indicating the assay to use to
#' obtain expression values. Must refer to a matrix-like object with integer or
#' numeric values. 
#' @param use.coldata Logical scalar indicating whether column metadata of x
#' should be included. Alternatively, a character or integer vector specifying
#' the column metadata fields to use. 
#' @param swap.rownames String or integer specifying the
#' [rowData][SummarizedExperiment::rowData] column containing the features
#' names. If `NULL`, `rownames(x)` will be used. 
#' @param check.names Logical scalar indicating whether column names of the
#' output should be made syntactically valid and unique. 
#' @rdname makePerCellDF
#' @export
methods::setMethod("makePerCellDF", "SummarizedExperiment", function(x, features = NULL, assay.type = NULL, use.coldata = TRUE, swap.rownames = NULL, check.names = FALSE) {
    # Initialize output list
    output <- .harvest_se_by_column(x,
        features = features,
        assay.type = assay.type,
        swap.rownames = swap.rownames
    )

    # Checking the names.
    output <- do.call(cbind, output)
    if (check.names) {
        colnames(output) <- make.names(colnames(output), unique = TRUE)
    }
    output
})

.harvest_se_by_column <- function(x, features, assay.type, swap.rownames = NULL) {
    # Collecting feature data from assay
    if (is.null(swap.rownames)) {
        all_feats <- rownames(x)
    } else {
        rowdata <- SummarizedExperiment::rowData(x)
        swap.rownames <- .use_names_to_integer_indices(
            swap.rownames, colnames(rowdata)
        )
        all_feats <- rowdata[[swap.rownames]]
    }
    features <- .use_names_to_integer_indices(features, all_feats)
    if (length(features)) {
        assay_data <- SummarizedExperiment::assay(x, assay.type)[
            features, ,
            drop = FALSE
        ]
        assay_data <- data.frame(as.matrix(t(assay_data)),
            check.names = FALSE, fix.empty.names = FALSE,
            row.names = colnames(x)
        )
        colnames(assay_data) <- all_feats[features]
    } else {
        assay_data <- NULL
    }

    # Collecting the column metadata.
    coldata <- SummarizedExperiment::colData(x)
    use.coldata <- .use_names_to_integer_indices(
        use.coldata,
        names = colnames(coldata)
    )
    if (length(use.coldata)) {
        cd <- coldata[, use.coldata, drop = FALSE]
        coldata <- data.frame(cd, check.names = FALSE, fix.empty.names = FALSE)
    } else {
        coldata <- NULL
    }

    list(assay = assay_data, coldata = coldata)
}
