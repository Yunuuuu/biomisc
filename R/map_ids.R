######################################################################
#' Gets the mapped ids (column) for a set of keys that are of a particular
#' keytype
#'
#' @description This function maps the feature names to a specific keytype.
#' Unmapped features are removed, and any duplicate mapped features are also
#' removed, retaining only the first occurrence based on the mean values across
#' all samples. 
#' @param x A matrix-like object.
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
#' @param keytype The keytype that matches the keys (`rownames(x)`) used.
#' @param decreasing A bool, Should the sort order be increasing or decreasing?
#' @rdname map_ids
#' @export
methods::setMethod("map_ids", "ANY", function(x, db, column, keytype, decreasing = TRUE) {
    feature_data <- map_ids_internal(
        id = rownames(x), db = db,
        column = column, keytype = keytype,
        decreasing = decreasing,
        order_by = rowMeans(x, na.rm = TRUE)
    )
    x <- x[feature_data$idx, , drop = FALSE]
    rownames(x) <- feature_data[[column]]
    x
})

#' @param rowdata_name Use column in `rowData(x)` to map feature ids, if `NULL`,
#' the `rownames(x)` will be used. 
#' @param assay_name A string or integer scalar indicating which `assays` in the
#' `x` to calculate the mean value for each row. IF `NULL`, the first assay will
#' be used.
#' @rdname map_ids
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
methods::setMethod("map_ids", "SummarizedExperiment", function(x, db, column, keytype, rowdata_name = NULL, assay_name = NULL, decreasing = TRUE) {
    assert_(rowdata_name, function(x) {
        length(x) == 1L && (is.character(x) || is.numeric(x))
    }, "a string or integer", null_ok = TRUE)
    if (is.null(rowdata_name)) {
        id <- rownames(x)
    } else {
        id <- SummarizedExperiment::rowData(x)[[rowdata_name]]
    }
    feature_data <- map_ids_internal(
        id = id, db = db,
        column = column, keytype = keytype, decreasing = decreasing,
        order_by = rowMeans(
            SummarizedExperiment::assay(x, i = assay_name %||% 1L),
            na.rm = TRUE
        )
    )
    x <- x[feature_data$idx, ]
    rownames(x) <- feature_data[[column]]
    SummarizedExperiment::rowData(x) <- cbind(
        feature_data, SummarizedExperiment::rowData(x)
    )
    x
})

map_ids_internal <- function(id, db, column, keytype, order_by, decreasing = TRUE, call = rlang::caller_env()) {
    assert_pkg("AnnotationDbi", fun = "map_ids", call = call)
    out <- data.table::substitute2(
        data.table::data.table(name = id), list(name = keytype)
    )
    out$idx <- seq_along(id)
    out[[column]] <- AnnotationDbi::mapIds(
        x = db, keys = out[[keytype]], column = column,
        keytype = keytype
    )
    out$val <- order_by
    expr <- data.table::substitute2(
        out[x != "" & !is.na(x)][order(val, decreasing = decreasing)], # nolint
        list(x = column)
    )
    out <- eval(expr)
    out <- data.table::unique(out, by = column, cols = c("idx", keytype))
    data.table::setcolorder(out, c("idx", keytype, column))
    data.table::setDF(out, rownames = out[[column]])
}
utils::globalVariables("val")
