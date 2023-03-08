`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

has_names <- function(x) {
    nms <- names(x)
    if (is.null(nms)) {
        return(rep_len(FALSE, length(x)))
    }
    !is.na(nms) & nms != ""
}

#' Report if an argument is a specific class
#'
#' @keywords internal
#' @noRd
assert_class <- function(x, is_class, class, null_ok = FALSE, arg = rlang::caller_arg(x), call = parent.frame()) {
    message <- "{.cls {class}} object"
    if (null_ok) {
        message <- paste(message, "or {.code NULL}", sep = " ")
    }
    message <- sprintf("{.arg {arg}} must be a %s", message)
    if (rlang::is_scalar_character(is_class)) {
        is_class <- substitute(function(x) {
            inherits(x, what = what)
        }, list(what = is_class))
        is_class <- eval(is_class)
    }

    # class sometimes can also contain the `NULL`
    if (is.null(x) && !is_class(x)) {
        if (!null_ok) {
            cli::cli_abort(c(message,
                "x" = "You've supplied a {.code NULL}"
            ), call = call)
        }
    } else if (!is_class(x)) {
        cli::cli_abort(c(message,
            "x" = "You've supplied a {.cls {class(x)}} object"
        ), call = call)
    }
}

#' Report if an argument has specific length
#' @keywords internal
#' @noRd
assert_length <- function(x, length, null_ok = FALSE, arg = rlang::caller_arg(x), call = parent.frame()) {
    length <- as.integer(length)
    if (length == 1L) {
        message <- "{.field scalar} object"
    } else {
        message <- "length {.val {length}} object"
    }
    if (null_ok) {
        message <- paste(message, "or {.code NULL}", sep = " ")
    }
    message <- sprintf("{.arg {arg}} must be a %s", message)
    is_right_length <- length(x) == length
    if (is.null(x) && !is_right_length) {
        if (!null_ok) {
            cli::cli_abort(c(message,
                "x" = "You've supplied a {.code NULL}"
            ), call = call)
        }
    } else if (!is_right_length) {
        cli::cli_abort(c(message,
            "x" = "You've supplied a length {.val {length(x)}} object"
        ), call = call)
    }
}

assert_pkg <- function(pkg, fun = NULL, call = parent.frame()) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (is.null(fun)) {
            call <- rlang::eval_bare(quote(match.call()), env = parent.frame())
            fun <- rlang::as_string(call[[1L]])
        }
        cli::cli_abort(
            "{.pkg pkg} must be installed to use {.fn {fun}}.",
            call = call
        )
    }
}

is_scalar_numeric <- function(x) {
    length(x) == 1L && is.numeric(x)
}
