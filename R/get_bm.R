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
#' @return A `data.frame`
get_bm <- function(mart, attributes, split = "chromosome_name", ..., times = 10L) {
    assert_pkg("biomaRt")
    if (!rlang::is_scalar_character(split)) {
        cli::cli_abort("{.arg split} must be a scalar string")
    }
    split_list <- get_bm_safe(
        mart = mart, attributes = split,
        ..., times = times
    )
    out <- lapply(cli::cli_progress_along(split_list), function(idx) {
        get_bm_safe(
            mart = mart,
            filters = split,
            attributes = attributes,
            values = split_list[[idx]],
            ...,
            times = times
        )
    })
    out <- data.table::rbindlist(out, use.names = FALSE)
    data.table::setDF(out)
    out
}

#' @keywords internal
#' @noRd
get_bm_safe <- function(..., times = 10L) {
    i <- 1L
    get_bm_rec <- function(...) {
        tryCatch(biomaRt::getBM(...), error = function(cnd) {
            if (i <= times) {
                i <<- i + 1L
                Recall(...)
            } else {
                cli::cli_abort("Cannot get data from {.field biomaRt} after trying {.val {times}} time{s?}")
            }
        })
    }
    get_bm_rec(...)
}
