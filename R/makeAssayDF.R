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

#' @param features Logical scalar indicating whether feature assay data of `x`
#' should be included. Alternatively, Character or integer vector specifying the
#' features for which to extract expression profiles across samples.
#' @param assay.type String or integer scalar indicating the assay to use to
#' obtain expression values. Must refer to a matrix-like object with integer or
#' numeric values. If `NULL`, the `1st` assay
#' (if x is a [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment]
#' object) or [exprs][Biobase::exprs] (if x is an
#' [ExpressionSet][Biobase::ExpressionSet] object) will be used.
#' @param use.coldata Logical scalar indicating whether column metadata of `x`
#' should be included. Alternatively, a character or integer vector specifying
#' the column metadata fields to use.
#' @param swap.rownames String or integer specifying the
#' [rowData][SummarizedExperiment::rowData] column containing the features
#' names. If `NULL`, `rownames(x)` will be used.
#' @param check.names Logical scalar indicating whether column names of the
#' output should be made syntactically valid and unique.
#' @name makePerCellDF
NULL

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname makePerCellDF
#' @export
methods::setMethod("makePerCellDF", "SummarizedExperiment", function(x, features = TRUE, assay.type = NULL, use.coldata = TRUE, swap.rownames = NULL, check.names = FALSE) {
    harvest_by_extra_data(x,
        features = features,
        use_extra_data = use.coldata,
        swap.rownames = swap.rownames,
        assay_data = SummarizedExperiment::assay(x, assay.type %||% 1L),
        rowdata = SummarizedExperiment::rowData(x),
        extra_data = SummarizedExperiment::colData(x),
        check.names = check.names,
        transpose = TRUE
    )
})

#' @importClassesFrom Biobase ExpressionSet
#' @rdname makePerCellDF
#' @export
methods::setMethod("makePerCellDF", "ExpressionSet", function(x, features = TRUE, assay.type = NULL, use.coldata = TRUE, swap.rownames = NULL, check.names = FALSE) {
    harvest_by_extra_data(x,
        features = features,
        use_extra_data = use.coldata,
        swap.rownames = swap.rownames,
        assay_data = Biobase::assayData(x)[[assay.type %||% "exprs"]],
        rowdata = Biobase::fData(x),
        extra_data = Biobase::pData(x),
        check.names = check.names,
        transpose = TRUE
    )
})

########################### makePerFeatureDF ###############################
#' Create a per-feature data.frame
#'
#' Create a per-feature data.frame (i.e., where each row represents a feature /
#' gene) from a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment], most
#' typically for creating custom ggplot2 plots.
#' @param x A [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment]
#'   or [ExpressionSet][Biobase::ExpressionSet] object.
#' @param ... Other arguments passed to specific methods.
#' @return A data.frame.
#' @name makePerFeatureDF
#' @export
methods::setGeneric(
    "makePerFeatureDF",
    function(x, ...) standardGeneric("makePerFeatureDF"),
    signature = "x"
)

#' @param use.rowdata Logical scalar indicating whether row metadata of `x`
#' should be included. Alternatively, a character or integer vector specifying
#' the row metadata fields to use.
#' @name makePerFeatureDF
NULL

#' @inheritParams makePerCellDF
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @rdname makePerFeatureDF
#' @export
methods::setMethod("makePerFeatureDF", "SummarizedExperiment", function(x, features = TRUE, assay.type = NULL, use.rowdata = TRUE, swap.rownames = NULL, check.names = FALSE) {
    rowdata <- SummarizedExperiment::rowData(x)
    harvest_by_extra_data(x,
        features = features,
        use_extra_data = use.rowdata,
        swap.rownames = swap.rownames,
        assay_data = SummarizedExperiment::assay(x, assay.type %||% 1L),
        rowdata = rowdata, extra_data = rowdata,
        check.names = check.names,
        transpose = FALSE
    )
})

#' @inheritParams makePerCellDF
#' @importClassesFrom Biobase ExpressionSet
#' @rdname makePerFeatureDF
#' @export
methods::setMethod("makePerFeatureDF", "ExpressionSet", function(x, features = TRUE, assay.type = NULL, use.rowdata = TRUE, swap.rownames = NULL, check.names = FALSE) {
    rowdata <- Biobase::fData(x)
    harvest_by_extra_data(x,
        features = features,
        use_extra_data = use.rowdata,
        swap.rownames = swap.rownames,
        assay_data = Biobase::assayData(x)[[assay.type %||% "exprs"]],
        rowdata = rowdata, extra_data = rowdata,
        check.names = check.names,
        transpose = FALSE
    )
})

########################### Utils Function ###############################
harvest_by_extra_data <- function(x, features, use_extra_data = TRUE, swap.rownames = NULL, check.names = FALSE, transpose, assay_data, rowdata, extra_data, call = rlang::caller_env()) {
    assay_data <- makeAssayDF(
        x = x, features = features,
        swap.rownames = swap.rownames, assay_data = assay_data,
        rowdata = rowdata,
        transpose = transpose,
        call = call
    )
    # Collecting the column metadata.
    use_extra_data <- use_names_to_integer_indices(
        use_extra_data,
        names = colnames(extra_data),
        call = call
    )
    if (length(use_extra_data)) {
        extra_data <- data.frame(extra_data,
            check.names = FALSE, fix.empty.names = FALSE
        )
        extra_data <- extra_data[, use_extra_data, drop = FALSE]
    } else {
        extra_data <- NULL
    }
    output <- cbind(assay_data, extra_data)
    # Checking the names.
    if (check.names) {
        colnames(output) <- make.names(colnames(output), unique = TRUE)
    }
    output
}

makeAssayDF <- function(x, features, swap.rownames, assay_data, rowdata, transpose, call = rlang::caller_env()) {
    assert_(swap.rownames, function(x) {
        length(x) == 1L && (is.character(x) || is.numeric(x))
    }, "a string or integer", null_ok = TRUE, call = call)

    # Collecting feature data from assay
    if (is.null(swap.rownames)) {
        all_feats <- rownames(x)
    } else {
        swap.rownames <- use_names_to_integer_indices(
            swap.rownames, colnames(rowdata),
            call = call
        )
        all_feats <- rowdata[[swap.rownames]]
    }
    features <- use_names_to_integer_indices(features, all_feats, call = call)
    if (length(features)) {
        assay_data <- assay_data[features, , drop = FALSE]
        if (transpose) assay_data <- t(assay_data)
        assay_data <- data.frame(as.matrix(assay_data),
            check.names = FALSE, fix.empty.names = FALSE
        )
        if (transpose) {
            rownames(assay_data) <- colnames(x)
            colnames(assay_data) <- all_feats[features]
        } else {
            rownames(assay_data) <- all_feats[features]
            colnames(assay_data) <- colnames(x)
        }
    } else {
        assay_data <- NULL
    }
    assay_data
}
