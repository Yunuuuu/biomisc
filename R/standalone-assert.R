# ---
# repo: Yunuuuu/biomisc
# file: standalone-assert.R
# last-updated: 2023-09-07
# license: https://unlicense.org
# dependencies: standalone-obj-type.R
# imports: rlang (>= 1.1.0)
# ---

#' Report if an argument is a specific class
#'
#' @param x The object type which does not conform to `what`. Its
#'   `obj_type_friendly()` is taken and mentioned in the error message.
#' @param what The friendly expected type as a string. Can be a
#'   character vector of expected types, in which case the error
#'   message mentions all of them in an "or" enumeration.
#' @param show_value Passed to `value` argument of `obj_type_friendly()`.
#' @param show_length Passed to `length` argument of `obj_type_friendly()`.
#' @param ... Arguments passed to [rlang::abort()].
#' @noRd
assert_ <- function(
    x, assert_fn, what,
    null_ok = FALSE,
    show_value = TRUE,
    show_length = FALSE,
    ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    if (null_ok && is.null(x) || assert_fn(x)) {
        return(invisible(NULL))
    }
    stop_input_type(x, what,
        null_ok = null_ok,
        show_value = show_value,
        show_length = show_length,
        ..., arg = arg, call = call
    )
}

stop_input_type <- function(
    x, what,
    null_ok = FALSE,
    show_value = TRUE,
    show_length = FALSE,
    ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    if (null_ok) {
        what <- c(what, format_code("NULL"))
    }
    if (length(what)) {
        what <- oxford_comma(what, final = "or")
    }
    msg <- sprintf(
        "%s must be %s, not %s.",
        format_arg(arg),
        what,
        obj_type_friendly(x, value = show_value, length = show_length)
    )
    rlang::abort(msg, ..., call = call, arg = arg)
}

assert_pkg <- function(pkg, version = NULL, fun = NULL, call = rlang::caller_env()) {
    if (!is_installed(pkg, version = version)) {
        if (is.null(fun)) {
            fun_call <- rlang::frame_call(frame = call)
            fun <- rlang::as_label(fun_call[[1L]])
        }
        pkg <- format_pkg(pkg)
        if (!is.null(version)) {
            pkg <- sprintf("%s (>=%s)", pkg, version)
        }
        rlang::abort(
            sprintf(
                "%s must be installed to use %s.",
                pkg, format_fn(fun)
            ),
            call = call
        )
    }
}

is_installed <- local({
    cache <- new.env()
    function(pkg, version = NULL) {
        id <- if (is.null(version)) pkg else paste(pkg, version, sep = ":")
        out <- cache[[id]]
        if (is.null(out)) {
            if (is.null(version)) {
                out <- requireNamespace(pkg, quietly = TRUE)
            } else {
                out <- requireNamespace(pkg, quietly = TRUE) &&
                    utils::packageVersion(pkg) >= version
            }
            cache[[id]] <<- out
        }
        out
    }
})

# scalar object
assert_string <- function(
    x, empty_ok = TRUE, na_ok = FALSE, show_length = TRUE, ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    what <- "a single string"
    if (!empty_ok) {
        what <- paste0(what, "(empty string is not allowed)")
    }
    if (na_ok) {
        what <- c(what, format_code("NA"))
    }
    assert_(
        x = x,
        assert_fn = function(x) {
            .rlang_check_is_string(x, empty_ok = empty_ok, na_ok = na_ok)
        }, what = what,
        show_length = show_length,
        ...,
        arg = arg,
        call = call
    )
}

.rlang_check_is_string <- function(x, empty_ok, na_ok) {
    if (rlang::is_string(x)) {
        if (empty_ok || !rlang::is_string(x, "")) {
            return(TRUE)
        }
    }
    if (na_ok && (identical(x, NA) || identical(x, NA_character_))) {
        return(TRUE)
    }
    FALSE
}

# S3 object
assert_s3_class <- function(
    x, is_class, what, ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    if (rlang::is_string(is_class)) {
        class <- is_class
        is_class <- function(x) {
            inherits(x, what = class)
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

assert_data_frame <- function(x, ..., arg = rlang::caller_arg(), call = rlang::caller_env()) {
    assert_s3_class(
        x = x, is_class = "data.frame",
        what = "a data frame",
        ..., arg = arg, call = call
    )
}

assert_data_frame_columns <- function(x, columns, ..., args = rlang::caller_arg(x), call = rlang::caller_env()) {
    missing_cols <- setdiff(columns, names(x))
    if (length(missing_cols)) {
        rlang::abort(
            sprintf(
                "One of %s must contain columns (%s), missing columns: %s",
                format_arg(args),
                oxford_comma(columns),
                oxford_comma(missing_cols)
            ), ...,
            call = call
        )
    }
}

# S4 object
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

stop_input_length <- function(
    x, length,
    scalar_ok = FALSE, null_ok = FALSE,
    ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    what <- NULL
    if (scalar_ok || length == 1L) {
        what <- c(what, sprintf("a %s", format_field("scalar")))
    }
    if (length != 1L) {
        what <- c(what, sprintf("of length %s", format_val(length)))
    }
    if (null_ok) {
        what <- c(what, format_code("NULL"))
    }
    if (length(what)) {
        what <- paste0(what, collapse = " or ")
    }
    msg <- sprintf(
        "%s must be %s, not of length %s.",
        format_arg(arg),
        what, format_val(length(x))
    )
    rlang::abort(msg, ..., call = call)
}


# Other utils
assert_data_frame_hierarchy <- function(x, parent_field, child_field = NULL, arg_children = rlang::caller_arg(child_field), ..., arg = rlang::caller_arg(x), call = rlang::caller_env()) {
    id1 <- format_code(sprintf("%s[[\"%s\"]]", arg, parent_field))
    id2 <- child_field %|n|%
        format_code(sprintf("%s[[\"%s\"]]", arg, child_field))
    assert_hierarchy(
        parents = x[[parent_field]],
        children = child_field %|n|% x[[child_field]],
        id1 = id1, id2 = id2, arg_children = arg_children,
        ..., call = call
    )
}

# assert hierarchy relationship, every child only has one parent
assert_hierarchy <- function(parents, children = NULL, id1 = rlang::caller_arg(parents), id2 = rlang::caller_arg(children), arg_children = rlang::caller_arg(children), ..., call = rlang::caller_env()) {
    if (is.null(children)) {
        n_unique <- length(unique(parents))
        if (n_unique > 1L) {
            rlang::abort(sprintf(
                "Only a unique value of %s can be used as no %s provided",
                id1, format_arg(arg_children, final = "or")
            ), ..., call = call)
        }
    } else {
        n_unique <- vapply(
            split(parents, children, drop = TRUE),
            function(x) length(unique(x)), integer(1L)
        )
        failed_idx <- n_unique > 1L
        if (any(failed_idx)) {
            rlang::abort(sprintf(
                "Multiple %s found in %s (%s)", id1, id2,
                format_val(names(n_unique)[failed_idx])
            ), ..., call = call)
        }
    }
}

assert_inclusive <- function(x, y, arg = rlang::caller_arg(x), call = rlang::caller_env()) {
    missing_items <- setdiff(x, y)
    if (length(missing_items)) {
        rlang::abort(
            sprintf(
                "Only values (%s) are allowed in %s, not %s",
                format_val(y), format_arg(arg), format_val(missing_items)
            ),
            call = call
        )
    }
}

#' The `format_` functions are easier to work with because they format the style
#' eagerly. However they produce slightly incorrect style in corner cases
#' because the formatting doesn't take into account the message type. In
#' principle, cli themes can create different stylings depending on the message
#' type.
#' @noRd 
format_val <- function(x, ...) .rlang_cli_format_inline(x, ".val", NULL, ...)
format_emph <- function(x, ...) .rlang_cli_format_inline(x, "emph", "%s", ...)
format_strong <- function(x, ...) {
    .rlang_cli_format_inline(x, "strong", "%s", ...)
}
format_code <- function(x, ...) .rlang_cli_format_inline(x, "code", "`%s`", ...)
format_q <- function(x, ...) .rlang_cli_format_inline(x, "q", NULL, ...)
format_pkg <- function(x, ...) .rlang_cli_format_inline(x, "pkg", NULL, ...)
format_fn <- function(x, ...) .rlang_cli_format_inline(x, "fn", "`%s()`", ...)
format_arg <- function(x, ...) .rlang_cli_format_inline(x, "arg", "`%s`", ...)
format_kbd <- function(x, ...) .rlang_cli_format_inline(x, "kbd", "[%s]", ...)
format_key <- function(x, ...) .rlang_cli_format_inline(x, "key", "[%s]", ...)
format_file <- function(x, ...) .rlang_cli_format_inline(x, "file", NULL, ...)
format_path <- function(x, ...) .rlang_cli_format_inline(x, "path", NULL, ...)
format_email <- function(x, ...) .rlang_cli_format_inline(x, "email", NULL, ...)
format_url <- function(x, ...) .rlang_cli_format_inline(x, "url", "<%s>", ...)
format_var <- function(x, ...) .rlang_cli_format_inline(x, "var", "`%s`", ...)
format_envvar <- function(x, ...) {
    .rlang_cli_format_inline(x, "envvar", "`%s`", ...)
}
format_field <- function(x, ...) .rlang_cli_format_inline(x, "field", NULL, ...)
format_cls <- function(x, ...) {
    fallback <- function(x) sprintf("<%s>", paste0(x, collapse = "/"))
    .rlang_cli_format_inline(x, "cls", fallback, ...)
}

#' @return A string
#' @noRd
.rlang_cli_format_inline <- function(x, span, fallback = "`%s`", ...) {
    if (inherits(x, "AsIs")) {
        oxford_comma(chr = x, ...)
    } else if (.rlang_cli_has_cli()) {
        cli::format_inline(oxford_comma(chr = paste0("{.", span, " {x}}"), ...))
    } else {
        oxford_comma(
            chr = .rlang_cli_style_inline(x, span, fallback = fallback), ...
        )
    }
}

#' @return An atomic character with the same length of x
#' @noRd
.rlang_cli_style_inline <- function(x, span, fallback = "`%s`") {
    if (.rlang_cli_has_cli()) {
        paste0("{.", span, " {\"", encodeString(x), "\"}}")
    } else if (is.null(fallback)) {
        x
    } else if (is.function(fallback)) {
        fallback(x)
    } else {
        sprintf(fallback, x)
    }
}

.rlang_cli_has_cli <- function(version = "3.0.0") {
    is_installed("cli", version = version)
}

oxford_comma <- function(chr, sep = ", ", final = "and") {
    n <- length(chr)

    if (n < 2L) {
        return(chr)
    }

    head <- chr[seq_len(n - 1L)]
    last <- chr[n]

    head <- paste(head, collapse = sep)

    # Write a or b. But a, b, or c.
    if (n > 2L) {
        paste0(head, sep, final, " ", last)
    } else {
        paste0(head, " ", final, " ", last)
    }
}
