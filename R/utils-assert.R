
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

assert_pkg <- function(pkg, fun = NULL, frame = rlang::caller_env()) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        if (is.null(fun)) {
            fun_call <- rlang::frame_call(frame = frame)
            fun <- rlang::as_label(fun_call[[1L]])
        }
        cli::cli_abort(
            "{.pkg {pkg}} must be installed to use {.fn {fun}}.",
            call = frame
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

# assert hierarchy relationship, every value only correspond to a unique id
assert_nest <- function(data, uid, group = NULL, arg = rlang::caller_arg(data), tag_uid = uid, tag_group = group, msg = NULL, cross_msg = NULL, cross_format = c("pair", "group", "uid"), info_msg = NULL, call = parent.frame()) {
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
