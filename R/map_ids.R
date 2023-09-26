######################################################################
#' Gets the mapped ids (column) for a set of keys that are of a particular
#' keytype
#'
#' @description This function maps the feature names to a specific keytype.
#' Unmapped features are removed, and any duplicate mapped features are also
#' removed, retaining only the first occurrence based on the mean values across
#' all samples.
#' @param x A matrix-like entity, such as a
#' [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment] or an
#' [ExpressionSet][Biobase::ExpressionSet] object. 
#' @param ... Other arguments passed to specific methods.
#' @name map_ids
#' @export
methods::setGeneric(
    "map_ids",
    function(x, ...) standardGeneric("map_ids"),
    signature = "x"
)

#' @param db An [AnnotationDb][AnnotationDbi::AnnotationDb] object. But in
#' practice this will mean an object derived from an AnnotationDb object such as
#' a `OrgDb` or `ChipDb` object.
#' @param column The columns or kinds of things that can be retrieved from the
#' database.
#' @param keytype The keytype that matches the keys (`rownames(x)` or column
#' specified in `swap_rownames`) used.
#' @param decreasing A bool, Should the sort order be increasing or decreasing?
#' @rdname map_ids
#' @export
methods::setMethod("map_ids", "ANY", function(x, db, column, keytype, decreasing = TRUE) {
    feature_data <- map_ids_internal(
        id = rownames(x), db = db,
        column = column, keytype = keytype,
        decreasing = decreasing,
        order_by_matrix = x
    )
    x <- x[feature_data$idx, , drop = FALSE]
    rownames(x) <- feature_data[[column]]
    x
})

#' @rdname map_ids
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
methods::setMethod("map_ids", "SummarizedExperiment", function(x, db, column, keytype, swap_rownames = NULL, assay_name = NULL, decreasing = TRUE) {
    out <- map_ids_swap(
        x = x, db = db, column = column, keytype = keytype,
        decreasing = decreasing, swap_rownames = swap_rownames,
        rowdata = SummarizedExperiment::rowData(x),
        assay = SummarizedExperiment::assay(x, i = assay_name %||% 1L)
    )
    SummarizedExperiment::rowData(out$x) <- cbind(
        out$feature_data, SummarizedExperiment::rowData(out$x)
    )
    out$x
})

#' @rdname map_ids
#' @importClassesFrom Biobase ExpressionSet
#' @export
methods::setMethod("map_ids", "ExpressionSet", function(x, db, column, keytype, swap_rownames = NULL, assay_name = NULL, decreasing = TRUE) {
    out <- map_ids_swap(
        x = x, db = db, column = column, keytype = keytype,
        decreasing = decreasing, swap_rownames = swap_rownames,
        rowdata = Biobase::fData(x),
        assay = Biobase::assayData(x)[[assay_name %||% "exprs"]]
    )
    Biobase::fData(out$x) <- cbind(out$feature_data, Biobase::fData(out$x))
    out$x
})

#' @param swap_rownames String or integer specifying the
#' [rowData][SummarizedExperiment::rowData] or [fData][Biobase::fData] column
#' containing the features names. If `NULL`, `rownames(x)` will be used.
#' @param assay_name A string or integer scalar indicating which assay in the
#' `x` to calculate the mean value for each row. If `NULL`, the `1st` assay (if
#' x is a [SummarizedExperiment][SummarizedExperiment::SummarizedExperiment]
#' object) or [exprs][Biobase::exprs] (if x is an
#' [ExpressionSet][Biobase::ExpressionSet] object) will be used.
#' @name map_ids
NULL

map_ids_swap <- function(x, db, column, keytype, decreasing, swap_rownames, rowdata, assay, call = rlang::caller_call()) {
    assert_(swap_rownames, function(x) {
        length(x) == 1L && (is.character(x) || is.numeric(x))
    }, "a string or integer", null_ok = TRUE, call = call)
    if (is.null(swap_rownames)) {
        id <- rownames(x)
    } else {
        swap_rownames <- use_names_to_integer_indices(
            swap_rownames, colnames(rowdata)
        )
        id <- rowdata[[swap_rownames]]
    }
    feature_data <- map_ids_internal(
        id = id, db = db,
        column = column, keytype = keytype, decreasing = decreasing,
        order_by_matrix = assay,
        call = call
    )
    x <- x[feature_data$idx, ]
    feature_data <- feature_data[c(keytype, column)]
    rownames(x) <- feature_data[[column]]
    list(x = x, feature_data = feature_data)
}

map_ids_internal <- function(id, db, column, keytype, order_by_matrix, decreasing = TRUE, call = rlang::caller_env()) {
    assert_pkg("AnnotationDbi", fun = "map_ids", call = call)
    assert_string(keytype, empty_ok = FALSE, call = call)
    assert_string(column, empty_ok = FALSE, call = call)
    out <- eval(data.table::substitute2(
        data.table::data.table(name = id), list(name = keytype)
    ))
    out$idx <- seq_along(id)
    out[[column]] <- AnnotationDbi::mapIds(
        x = db, keys = out[[keytype]], column = column,
        keytype = keytype
    )
    out$val <- rowMeans(order_by_matrix, na.rm = TRUE)
    out <- eval(data.table::substitute2(
        out[.id != "" & !is.na(.id)][ # nolint
            order(val, decreasing = decreasing)
        ],
        list(.id = column)
    ))
    out <- unique(out, by = column, cols = c("idx", keytype))
    data.table::setcolorder(out, c("idx", keytype, column))
    data.table::setDF(out, rownames = out[[column]])
}
utils::globalVariables(c("val", ".id"))
