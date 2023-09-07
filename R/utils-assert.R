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

assert_s3_class <- function(
    x, is_class, what, ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    if (rlang::is_scalar_character(is_class)) {
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

assert_s4_class <- function(
    x, is_class, what, ...,
    arg = rlang::caller_arg(x),
    call = rlang::caller_env()) {
    if (rlang::is_scalar_character(is_class)) {
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
    if (inherits(arg, "AsIs")) {
        .format_arg <- identity
    } else {
        .format_arg <- format_arg
    }
    msg <- sprintf(
        "%s must be %s, not %s.",
        .format_arg(arg),
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
    if (inherits(arg, "AsIs")) {
        .format_arg <- identity
    } else {
        .format_arg <- format_arg
    }
    msg <- sprintf(
        "%s must be %s, not of length %s.",
        .format_arg(arg),
        what, format_val(length(x))
    )
    rlang::abort(msg, ..., call = call)
}

assert_data_frame <- function(x, ..., arg = rlang::caller_arg(), call = rlang::caller_env()) {
    assert_s3_class(x, "data.frame",
        what = "a data frame",
        ..., arg = arg, call = call
    )
}

assert_data_frame_columns <- function(x, columns, ..., args = rlang::caller_arg(x), call = rlang::caller_env()) {
    missing_cols <- setdiff(columns, names(x))
    if (length(missing_cols)) {
        if (inherits(args, "AsIs")) {
            .format_arg <- oxford_comma
        } else {
            .format_arg <- format_arg
        }
        rlang::abort(
            sprintf(
                "Cannot find columns (%s) in %s",
                oxford_comma(missing_cols),
                .format_arg(args, final = "or")
            ), ...,
            call = call
        )
    }
}

# assert hierarchy relationship, every value only correspond to a unique id
assert_nest <- function(data, uid, group = NULL, arg = rlang::caller_arg(data), tag_uid = uid, tag_group = group, msg = NULL, cross_msg = NULL, cross_format = c("pair", "group", "uid"), info_msg = NULL, call = rlang::caller_env()) {
    if (is.null(group)) {
        pairs <- list(unique(uid))
    } else {
        pairs <- lapply(
            split(data[[uid]], data[[group]], drop = TRUE),
            unique
        )
    }
    pairs <- pairs[lengths(pairs) > 1L]
    if (length(pairs) > 0L) {
        cross_format <- match.arg(cross_format)
        if (is.null(group)) {
            msg <- msg %||% "A single {.field {tag_group}}"
            cross_msg <- cross_msg %||%
                switch(cross_format,
                    group = NULL,
                    uid = ,
                    pair = "Incorrect {.field {tag_uid}}: {.val {pairs[[1L]]}}",
                )
        } else {
            msg <- msg %||% "Every {.field {tag_group}}"
            cross_msg <- cross_msg %||%
                switch(cross_format,
                    group = "{tag_group} with multiple {tag_uid}: {.val {names(pairs)}}",
                    uid = "{tag_uid} with multiple {tag_group}: {.val {pairs}}",
                    pair = vapply(seq_along(pairs), function(i) {
                        sprintf(
                            "{.field {names(pairs)[[%d]]}}: {.val {pairs[[%d]]}}",
                            i, i
                        )
                    }, character(1L), USE.NAMES = FALSE)
                )
        }
        msg <- paste(
            msg,
            "can only have one {.field {tag_uid}} value in {.arg {arg}}"
        )
        names(cross_msg) <- rlang::rep_along(cross_msg, "x")
        cli::cli_abort(c(msg, cross_msg, i = info_msg), call = call)
    }
}

assert_in <- function(x, y, arg_x = rlang::caller_arg(x), call = rlang::caller_env()) {
    missing_items <- setdiff(x, y)
    if (length(missing_items)) {
        cli::cli_abort(
            c(
                "value allowed in {.arg {arg_x}}: {.val {y}}",
                x = "erroneous value{?s}: {.val {missing_items}}"
            ),
            call = call
        )
    }
}

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
    if (.rlang_cli_has_cli()) {
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
