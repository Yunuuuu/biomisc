#' Get AnnotationHub object
#' @param ... See [AnnotationHub][AnnotationHub::AnnotationHub] for additional
#' arguments.
#' @param version A valid Bioconductor version, see:
#' <https://bioconductor.org/config.yaml> for all version.
#' [AnnotationHub][AnnotationHub::AnnotationHub()] restrict the AnnotationHub
#' snapshotDate before the Bioconductor version, this function just bypass the
#' restriction to let you query all AnnotationHub data.
#' @param hub The URL for the online AnnotationHub.
#' @param cache The file system location of the local AnnotationHub cache.
#' @param proxy Set the proxy.
#' @param localHub A bool, whether to use the cache.
#' @return A [AnnotationHub][AnnotationHub::AnnotationHub] instance.
#' @export
anno_hub <- function(..., version = NULL, hub = NULL, cache = NULL, proxy = NULL, localHub = NULL) {
    assert_pkg("AnnotationHub")
    assert_pkg("BiocManager")
    version <- version %||% BiocManager::version()
    url <- hub %||% AnnotationHub::getAnnotationHubOption("URL")
    cache <- cache %||% AnnotationHub::getAnnotationHubOption("CACHE")
    proxy <- proxy %||% AnnotationHub::getAnnotationHubOption("PROXY")
    localHub <- localHub %||% AnnotationHub::getAnnotationHubOption("LOCAL")
    db_path <- AnnotationHub:::.create_cache(
        "AnnotationHub",
        url = url, cache, proxy, localHub,
        ask = AnnotationHub::getAnnotationHubOption("ASK")
    )
    if (!localHub) {
        tryCatch(
            {
                dates <- as.POSIXlt(
                    AnnotationHub:::.possibleDates(db_path),
                    format = "%Y-%m-%d"
                )
                restrict <- as.POSIXlt(
                    AnnotationHub:::.biocVersionDate(version),
                    format = "%Y-%m-%d"
                )
                if (length(restrict)) {
                    db_date <- as.character(max(dates[dates <= restrict]))
                } else {
                    db_date <- as.character(max(dates))
                }
            },
            error = function(err) {
                stop(
                    "failed to connect", "\n  reason: ", conditionMessage(err),
                    "\n  Consider rerunning with 'localHub=TRUE'"
                )
            }
        )
    } else {
        dates <- AnnotationHub:::.possibleDates(db_path)
        db_date <- dates[length(dates)]
    }
    db_uid <- AnnotationHub:::.db_uid0(db_path, db_date, localHub)
    hub <- methods::new("AnnotationHub",
        cache = cache, hub = url, date = db_date,
        .db_path = db_path, .db_uid = db_uid, isLocalHub = localHub,
        ...
    )
    message("snapshotDate(): ", AnnotationHub::snapshotDate(hub))
    if (!localHub) {
        index <- AnnotationHub:::.db_create_index(hub)
        AnnotationHub:::.db_index(hub) <- index
    } else {
        index <- AnnotationHub:::.db_index_file(hub)
        AnnotationHub:::.db_index(hub) <- index
        hub <- AnnotationHub:::.subsethub(hub)
    }
    hub
}
