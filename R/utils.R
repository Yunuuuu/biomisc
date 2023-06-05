`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

#' @importFrom data.table %chin%
NULL

#' Report if an argument is a specific class
#'
#' @keywords internal
#' @noRd
assert_class <- function(
    x, is_class, msg,
    cross_msg = "You've supplied a {.cls {class(x)}} object",
    null_ok = FALSE, arg = rlang::caller_arg(x), call = parent.frame(),
    .envir = environment()) {
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
            ), call = call, .envir = .envir)
        }
    } else if (!is_right_class) {
        if (!is.null(cross_msg)) {
            msg <- c(msg, x = cross_msg)
        }
        cli::cli_abort(msg, call = call, .envir = .envir)
    }
}

#' Report if an argument has specific length
#' @keywords internal
#' @noRd
assert_length <- function(x, length, msg, scalar_ok = FALSE, null_ok = FALSE, arg = rlang::caller_arg(x), call = parent.frame(), .envir = environment()) {
    if (!missing(length)) {
        length <- as.integer(length)
        if (missing(msg)) {
            if (length > 1L) {
                msg <- "of length {.val {length}}"
            } else if (length == 1L) {
                msg <- "a {.field scalar}"
            }
        }
    }
    if (scalar_ok) {
        if (missing(msg)) {
            msg <- "a {.field scalar}"
        } else if (!missing(length) && length > 1L) {
            msg <- paste("a {.field scalar} or", msg, sep = " ")
        }
    }
    if (null_ok) {
        msg <- paste(msg, "or {.code NULL}", sep = " ")
    }
    msg <- sprintf("{.arg {arg}} must be %s", msg)
    if (!missing(length)) {
        is_right_length <- length(x) == length
    } else {
        is_right_length <- FALSE
    }
    if (scalar_ok) {
        is_right_length <- is_right_length || length(x) == 1L
    }

    if (is.null(x) && !is_right_length) {
        if (!null_ok) {
            cli::cli_abort(c(msg,
                "x" = "You've supplied a {.code NULL}"
            ), call = call, .envir = .envir)
        }
    } else if (!is_right_length) {
        cli::cli_abort(c(msg,
            "x" = "You've supplied a length {.val {length(x)}} object"
        ), call = call, .envir = .envir)
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

assert_df_with_columns <- function(x, cols, check_class = length(arg) == 1L, arg = rlang::caller_arg(x), call = parent.frame()) {
    if (check_class) {
        is_right <- inherits(x, "data.frame")
        msg <- c(i = "{.arg {arg}} must be a {.cls data.frame}")
    } else {
        is_right <- TRUE
        msg <- character()
    }
    arg_style <- cli::cli_vec(arg, style = list( # nolint
        "vec-last" = " or ", "vec-trunc" = 10L
    ))
    cols_style <- cli::cli_vec(cols, style = list("vec-trunc" = 10L)) # nolint
    msg <- c(msg,
        i = "{.arg {arg_style}} must contatin {.val {cols_style}} column{?s}"
    )
    if (!is_right) {
        msg <- c(msg, x = "You've supplied a {.cls {class(x)}} object")
    }
    missing_cols <- setdiff(cols, names(x))
    is_right <- length(missing_cols) == 0L
    if (!is_right) {
        msg <- c(msg, x = "Cannot find column{?s}: {.val {missing_cols}}")
    }
    if (!is_right) cli::cli_abort(msg, call = call)
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

trim_value <- function(x, threshold = 1 - .Machine$double.neg.eps) {
    x[x >= threshold] <- threshold
    x
}
