KEGG_ORGANISM <- new.env()
KEGG_ORGANISM$has_data <- FALSE
kegg_organism <- function() {
    assert_pkg("KEGGREST")
    if (is.null(out <- .subset2(KEGG_ORGANISM, "has_data"))) {
        organism <- apply(
            KEGGREST::keggList("organism"),
            2L, identity,
            simplify = FALSE
        )
        list2env(organism, envir = KEGG_ORGANISM)
        assign("has_data", TRUE, envir = KEGG_ORGANISM)
    }
    out
}

#' Get database from KEGG
#'
#' @inheritParams KEGGREST::keggList
#' @return A [data.table][data.table::data.table].
#' @export
keggdb <- function(database, organism = "hsa") {
    assert_pkg("KEGGREST")
    database <- match.arg(database, KEGGREST::listDatabases())
    items <- KEGGREST::keggList(database, organism)
    out <- vector("list", length(items))
    identifiers <- names(items)
    names(out) <- identifiers
    index <- seq_along(out)
    groups <- split(index, (index %/% 10) + 1L)
    progress_id <- cli::cli_progress_bar(
        name = "keggGet",
        format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} [{cli::pb_rate}] | {cli::pb_eta_str}",
        format_done = "Get from KEGG {database} for {.val {cli::pb_total}} quer{?y/ies} in {cli::pb_elapsed}",
        total = length(out),
        clear = FALSE
    )
    for (idx in groups) {
        out[idx] <- KEGGREST::keggGet(identifiers[idx])
        cli::cli_progress_update(inc = length(idx), id = progress_id)
    }
    out <- lapply(out, function(kegg_list) {
        data.table::as.data.table(lapply(kegg_list, function(item) {
            if (length(item) > 1L) list(item) else item
        }))
    })
    data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

kegg_query <- function(..., database, option = NULL) {
    database <- match.arg(database, KEGGREST::listDatabases())
    if (is.null(option)) {
        KEGGREST::keggFind(database = database, query = c(...))
    } else {
        KEGGREST::keggFind(database = database, query = c(...), option = option)
    }
}
