`%||%` <- function(x, y) {
    if (is.null(x)) y else x
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

read_internal_extdata <- function(...) {
    readRDS(system.file("extdata", ..., package = "biomisc"))
}

trim_value <- function(x, threshold = 1 - .Machine$double.neg.eps) {
    x[x >= threshold] <- threshold
    x
}
