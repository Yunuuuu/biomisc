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
#' a `OrgDb` or `ChipDb` object. Alternatively, a named atomic vector can be
#' used to match `keytype` by names.
#' @param column The columns or kinds of things that can be retrieved from the
#' `db` database.
#' @param keytype The keytype that matches the keys (`rownames(x)` or column
#' specified in `swap_rownames`) used.
#' @param decreasing A bool, Should the sort order be increasing or decreasing?
#' @rdname map_ids
#' @export
methods::setMethod("map_ids", "ANY", function(x, db, column, keytype, decreasing = TRUE) {
    out <- map_ids_internal(
        x = x, db = db,
        column = column, keytype = keytype,
        decreasing = decreasing,
        assay = x
    )
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

#' @rdname map_ids
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
methods::setMethod("map_ids", "SummarizedExperiment", function(x, db, column, keytype, swap_rownames = NULL, assay_name = NULL, decreasing = TRUE) {
    out <- map_ids_internal(
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
    out <- map_ids_internal(
        x = x, db = db, column = column, keytype = keytype,
        decreasing = decreasing, swap_rownames = swap_rownames,
        rowdata = Biobase::fData(x),
        assay = Biobase::assayData(x)[[assay_name %||% "exprs"]]
    )
    Biobase::fData(out$x) <- cbind(out$feature_data, Biobase::fData(out$x))
    out$x
})

map_ids_internal <- function(x, db, column, keytype, decreasing, assay, swap_rownames = NULL, rowdata, call = rlang::caller_env()) {
    assert_pkg("AnnotationDbi", fun = "map_ids", call = call)
    assert_string(keytype, empty_ok = FALSE, call = call)
    assert_string(column, empty_ok = FALSE, call = call)
    assert_bool(decreasing, call = call)
    assert_(swap_rownames, function(x) {
        length(x) == 1L && (is.character(x) || is.numeric(x))
    }, "a string or integer", null_ok = TRUE, call = call)

    # swap rownames ------------------------
    if (is.null(swap_rownames)) {
        id <- rownames(x)
        if (is.null(id)) {
            rlang::abort(sprintf("no rownames found in %s", style_arg(x)),
                call = call
            )
        }
    } else {
        swap_rownames <- use_names_to_integer_indices(
            swap_rownames, colnames(rowdata),
            call = call
        )
        id <- rowdata[[swap_rownames]]
    }

    # map ids -------------------------------
    feature_data <- eval(data.table::substitute2(
        data.table::data.table(name = id), list(name = keytype)
    ))
    feature_data$idx <- seq_along(id)
    if (rlang::is_atomic(db)) {
        if (rlang::is_named(db)) {
            feature_data[[column]] <- db[feature_data[[keytype]]]
        } else {
            rlang::abort(
                sprintf("atomic %s must be named", style_arg("db")),
                call = call
            )
        }
    } else {
        feature_data[[column]] <- AnnotationDbi::mapIds(
            x = db, keys = feature_data[[keytype]], column = column,
            keytype = keytype
        )
    }

    # filter duplicates ------------------------------
    feature_data$.val <- rowMeans(assay, na.rm = TRUE)
    feature_data <- eval(data.table::substitute2(
        feature_data[.id != "" & !is.na(.id)], # nolint,
        list(.id = column)
    ))
    if (decreasing) {
        data.table::setorderv(feature_data, c(column, ".val"), c(1L, -1L))
    } else {
        data.table::setorderv(feature_data, c(column, ".val"), c(1L, 1L))
    }
    feature_data <- unique(feature_data, by = column)
    x <- x[feature_data$idx, , drop = FALSE]
    feature_data <- feature_data[, .SD, .SDcols = c(keytype, column)]
    data.table::setDF(feature_data, rownames = feature_data[[column]])
    list(x = x, feature_data = feature_data)
}

utils::globalVariables(".id")
