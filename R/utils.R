`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

#' @importFrom data.table %chin%
NULL

#' Report if an argument is a specific class
#'
#' @keywords internal
#' @noRd
assert_class <- function(x, is_class, msg, cross_msg = "{.cls {class(x)}} object", null_ok = FALSE, arg = rlang::caller_arg(x), call = parent.frame()) {
    if (rlang::is_scalar_character(is_class)) {
        class <- is_class
        is_class <- function(x) {
            inherits(x, what = class)
        }
        if (missing(msg)) {
            msg <- "{.cls {class}} object"
        }
    }
    if (null_ok) {
        msg <- paste(msg, "or {.code NULL}", sep = " ")
    }
    msg <- sprintf("{.arg {arg}} must be a %s", msg)
    is_right_class <- is_class(x)
    # is_class sometimes return `TRUE` for`NULL`
    if (is.null(x) && !is_right_class) {
        if (!null_ok) {
            cli::cli_abort(c(msg,
                "x" = "You've supplied a {.code NULL}"
            ), call = call)
        }
    } else if (!is_right_class) {
        if (!is.null(cross_msg)) {
            cross_msg <- sprintf("You've supplied a ", cross_msg)
            msg <- c(msg, x = cross_msg)
        }
        cli::cli_abort(msg, call = call)
    }
}

#' Report if an argument has specific length
#' @keywords internal
#' @noRd
assert_length <- function(x, length, msg, scalar_ok = FALSE, null_ok = FALSE, arg = rlang::caller_arg(x), call = parent.frame()) {
    if (missing(length)) {
        length <- as.integer(length)
        if (missing(msg)) {
            if (length > 1L) {
                msg <- "length {.field {length}}"
            } else if (length == 1L) {
                msg <- "a {.field scalar}"
            }
        }
    }
    if (scalar_ok) {
        if (missing(msg)) {
            msg <- "a {.field scalar}"
        } else if (!missing(length) && length > 1L) {
            msg <- paste(msg, "or a {.field scalar}", sep = " ")
        }
    }
    if (null_ok) {
        msg <- paste(msg, "or {.code NULL}", sep = " ")
    }
    msg <- sprintf("{.arg {arg}} must be of %s", msg)
    is_right_length <- length(x) == length
    if (is.null(x) && !is_right_length) {
        if (!null_ok) {
            cli::cli_abort(c(msg,
                "x" = "You've supplied a {.code NULL}"
            ), call = call)
        }
    } else if (!is_right_length) {
        cli::cli_abort(c(msg,
            "x" = "You've supplied a length {.val {length(x)}} object"
        ), call = call)
    }
}

assert_pkg <- function(pkg, fun = NULL, call = parent.frame()) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (is.null(fun)) {
            fn_call <- eval(quote(match.call()), envir = parent.frame())
            fun <- rlang::as_string(fn_call[[1L]])
        }
        cli::cli_abort(
            "{.pkg {pkg}} must be installed to use {.fn {fun}}.",
            call = call
        )
    }
}

is_scalar_numeric <- function(x) {
    length(x) == 1L && is.numeric(x)
}

format_percent <- function(x, digits = 2L, format = "f", ...) {
    sprintf("%s%%", formatC(x * 100L, format = format, digits = digits, ...))
}

read_internal_extdata <- function(...) {
    readRDS(system.file("extdata", ..., package = "biomisc"))
}
