assert_data_frame_columns <- function(x, columns, ..., args = rlang::caller_arg(x), call = rlang::caller_env()) {
    missing_cols <- setdiff(columns, names(x))
    if (length(missing_cols)) {
        args <- style_arg(args)
        if (length(args) == 1L) {
            msg <- args
        } else {
            msg <- sprintf("One of %s", oxford_comma(args, final = "or"))
        }
        rlang::abort(
            c(
                sprintf(
                    "%s must contain columns: %s", msg,
                    oxford_comma(columns)
                ),
                x = sprintf("missing columns: %s", oxford_comma(missing_cols))
            ),
            ...,
            call = call
        )
    }
}

# S4 object -------------------------------------
assert_s4_class <- function(
    x, is_class, what, ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    if (rlang::is_string(is_class)) {
        class <- is_class
        is_class <- function(x) {
            methods::is(x, class)
        }
        if (missing(what)) {
            what <- class
        }
    }
    assert_(
        x = x, assert_fn = is_class, what = what,
        ..., arg = arg, call = call
    )
}

#' Report if an argument has specific length
#' @keywords internal
#' @noRd
assert_length <- function(x, length, scalar_ok = FALSE, null_ok = FALSE, ..., arg = rlang::caller_arg(x), call = rlang::caller_env()) {
    is_right_length <- length(x) == length
    if (scalar_ok) {
        is_right_length <- is_right_length || length(x) == 1L
    }
    if (null_ok) {
        is_right_length <- is_right_length || is.null(x)
    }
    if (!is_right_length) {
        stop_input_length(x, length,
            scalar_ok = scalar_ok, null_ok = null_ok,
            ..., arg = arg, call = call
        )
    }
}

# Other assert function ---------------------------------
assert_data_frame_hierarchy <- function(x, parent_field, child_field = NULL, arg_children = rlang::caller_arg(child_field), ..., arg = rlang::caller_arg(x), call = rlang::caller_env()) {
    id_parent <- style_code(sprintf("%s[[\"%s\"]]", arg, parent_field))
    if (is.null(child_field)) {
        id_child <- NULL
    } else {
        id_child <- style_code(sprintf("%s[[\"%s\"]]", arg, child_field))
    }
    if (is.null(child_field)) {
        children <- NULL
    } else {
        children <- x[[child_field]]
    }
    assert_hierarchy(
        parents = x[[parent_field]],
        children = children,
        id_parent = id_parent, id_child = id_child, arg_children = arg_children,
        ..., call = call
    )
}

# assert hierarchy relationship, every child only has one parent
assert_hierarchy <- function(parents, children = NULL, id_parent = rlang::caller_arg(parents), id_child = rlang::caller_arg(children), arg_children = rlang::caller_arg(children), ..., call = rlang::caller_env()) {
    if (is.null(children)) {
        n_unique <- length(unique(parents))
        if (n_unique > 1L) {
            rlang::abort(sprintf(
                "Only a unique value of %s can be used as no %s provided",
                id_parent, oxford_comma(style_arg(arg_children), final = "or")
            ), ..., call = call)
        }
    } else {
        n_unique <- vapply(
            split(parents, children, drop = TRUE),
            function(x) length(unique(x)), integer(1L)
        )
        failed_idx <- n_unique > 1L
        if (any(failed_idx)) {
            rlang::abort(c(
                sprintf("Multiple %s found in %s", id_parent, id_child),
                x = sprintf(
                    "wrong %s: %s", id_child,
                    oxford_comma(style_val(names(n_unique)[failed_idx]))
                )
            ), ..., call = call)
        }
    }
}
