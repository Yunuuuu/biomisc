#' Retrieves information from the BioMart database
#' @description
#'   A modified [getBM][biomaRt::getBM], which will split our task into pieces
#'   to overcome the network restrictions.
#' @inheritParams biomaRt::getBM
#' @param split The attribute to split the task, since `get_bm` use filters to
#'   implement task spliting, user should ensure the split is in the filter
#'   attributes. You can list all filter attributes for a mart by
#'   [listFilters][biomaRt::listFilters].
#' @param ... Other arguments passed to [getBM][biomaRt::getBM].
#' @param times The maximal times to try to download data from biomaRt database.
#' @seealso [getBM][biomaRt::getBM]
#' @return A [data.table][data.table::data.table]
#' @export
get_bm <- function(mart, attributes, split = "chromosome_name", ..., times = 10L) {
    assert_pkg("biomaRt")
    assert_pkg("curl")
    if (!rlang::is_scalar_character(split)) {
        cli::cli_abort("{.arg split} must be a scalar string")
    }
    dots <- list(...)
    dots$curl <- NULL
    if (!is.null(dots$values)) {
        cli::cli_abort(c(
            "{.arg values} argument cannot work",
            "Use a named {.arg filters} instead"
        ))
    }
    split_list <- do.call(
        "get_bm_safe",
        c(dots, list(
            mart = mart, attributes = split,
            times = times
        ))
    )[[split]]
    out <- lapply(
        cli::cli_progress_along(split_list, "Downloading"),
        function(idx) {
            args <- c(dots, list(
                mart = mart, attributes = attributes,
                times = times, curl = curl::new_handle()
            ))
            args$filters <- c(
                args$filters,
                rlang::inject(rlang::list2(!!split := split_list[[idx]]))
            )
            do.call("get_bm_safe", args)
        }
    )
    data.table::rbindlist(out, use.names = FALSE)
}

#' @keywords internal
#' @noRd
get_bm_safe <- function(..., times = 10L) {
    i <- 1L
    get_bm_rec <- function(...) {
        tryCatch(biomaRt::getBM(...), error = function(cnd) {
            if (i <= times) {
                i <<- i + 1L
                get_bm_rec(...)
            } else {
                cli::cli_abort("Cannot get data from {.field biomaRt} after trying {.val {times}} time{?s}")
            }
        })
    }
    get_bm_rec(...)
}
