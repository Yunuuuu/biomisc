`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

`%|n|%` <- function(x, y) {
    if (is.null(x)) x else y
}

#' @importFrom data.table %chin%
NULL

silent_expr <- function(expr, msg, ...) {
    withCallingHandlers(
        rlang::eval_tidy(rlang::enquo(expr)),
        warning = function(cnd) {
            if (grepl(pattern = msg, x = conditionMessage(cnd), ...)) {
                rlang::cnd_muffle(cnd)
            }
        }
    )
}

is_scalar_numeric <- function(x) {
    length(x) == 1L && is.numeric(x)
}

format_percent <- function(x, digits = 2L, format = "f", ...) {
    sprintf("%s%%", formatC(x * 100L, format = format, digits = digits, ...))
}

unique_n <- function(x) {
    length(unique(x))
}

read_internal_extdata <- function(...) {
    readRDS(system.file("extdata", ..., package = "biomisc"))
}

use_names_to_integer_indices <- function(use, names, arg = rlang::caller_arg(use), call = rlang::caller_env()) {
    if (anyNA(use)) {
        rlang::abort(
            sprintf("%s cannot contain `NA`", style_arg(arg)),
            call = call
        )
    }
    if (isTRUE(use)) {
        use <- seq_along(names)
    } else if (isFALSE(use)) {
        use <- integer(0L)
    } else if (is.character(use)) {
        use <- match(use, names)
        if (anyNA(use)) {
            rlang::abort(sprintf(
                "%s contains invalid values", style_arg(arg)
            ), call = call)
        }
    } else if (is.numeric(use)) {
        if (any(use < 1L) || any(use > length(names))) {
            rlang::abort(sprintf(
                "%s contains out-of-bounds indices", style_arg(arg)
            ), call = call)
        }
    } else {
        rlang::abort(
            sprintf("%s must be a bool or an atomic numeic/character", style_arg(arg)),
            call = call
        )
    }
    use
}

trim_value <- function(x, threshold = 1 - .Machine$double.neg.eps) {
    x[x > threshold] <- threshold
    x
}

counts_matrix <- function(x, levels = NULL, y = NULL) {
    if (!is.null(levels)) {
        x <- factor(x, levels = levels)
    }
    if (is.null(y)) {
        counts_mat <- c(table(x))
        counts_mat <- matrix(counts_mat, ncol = 1L)
        rownames(counts_mat) <- levels(x)
    } else {
        counts_mat <- table(x, y)
        # convert table into matrix
        names(dimnames(counts_mat)) <- NULL
        counts_mat <- unclass(counts_mat)
    }
    counts_mat
}

cli_list <- function(list, ...) {
    main_item_lid <- cli::cli_ul()
    .mapply(function(name, value) {
        cli_named_li(name = name, value = value, ...)
    }, list(name = names(list), value = list), MoreArgs = NULL)
    cli::cli_end(id = main_item_lid)
}

cli_named_li <- function(
    name, value, label = "", sep = ": ", add_num = TRUE,
    cli_style = list("vec-trunc" = 3L)) {
    cli_value <- cli::cli_vec(value, cli_style) # nolint
    msg <- "{.field {name}}{sep}{.val {cli_value}}"
    if (add_num) {
        if (nzchar(label)) {
            suffix <- "({length(values)} {label}{?s})"
        } else {
            suffix <- "({length(values)})"
        }
    } else {
        if (nzchar(label)) {
            suffix <- "{label}{?s}"
        } else {
            suffix <- NULL
        }
    }
    msg <- paste(msg, suffix, sep = " ")
    cli::cli_li(msg)
}
