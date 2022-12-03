#' Seach CellMarker database
#' Details see http://xteam.xbio.top/CellMarker/index.jsp
#' @param markers the markers to search in the CellMarker database, should be
#' the gene symbol.
#' @param species a scalar string, "human" or "mouse".
#' @param internal logical value, indicates whether to use internal CellMarker
#' data or not. `NULL` indicates this will be determined by the CellMarker has
#' been downloaded once, if we have used this function once with a `internal`
#' FALSE, the value `NULL` will indicate `TRUE`. The internal data is downloaded
#' from CellMarker (2022-12-04).
#' @return a data.frame containing all search results, a column named `targeted`
#' including the intersection between `markers` (provided by the user) and the
#' cellMarker or geneSymbol column in CellMarker data.
#' @export
cellmarker_search <- function(markers, species = "human", internal = NULL) {
    data <- data.table::copy(cellmarker_get(species, internal))
    data[, targeted := lapply(gene_list, intersect, markers)] # nolint
    data.table::setcolorder(data, "targeted", before = "cellMarker")
    geneid_cols <- intersect(
        c("cellMarker", "geneSymbol", "geneID", "proteinName", "proteinID"),
        names(data)
    )
    data.table::setcolorder(
        data, geneid_cols,
        after = "gene_list"
    )
    data <- data[vapply(targeted, function(x) length(x) > 0L, logical(1L))] # nolint
    data.table::setDF(data)[]
}

cellmarker_get <- function(species = "human", internal = NULL) {
    species <- match.arg(species, c("human", "mouse"))
    if (is.null(internal)) {
        if (is.null(cellmarker_database_external[[species]])) {
            internal <- TRUE
        } else {
            internal <- FALSE
        }
    }
    if (internal) {
        cellmarker_database[[species]] # nolint
    } else {
        cellmarker_download(species)
    }
}

cellmarker_prepare <- function(data) {
    data[, gene_list := .mapply( # nolint
        union, unname(lapply(.SD, function(markers) {
            markers_list <- strsplit(
                gsub("\\s*\\[\\s*|\\s*\\]\\s*", "", markers, perl = TRUE),
                "\\s*,\\s*",
                perl = TRUE
            )
            lapply(markers_list, function(markers_trim) {
                markers_trim[
                    !is.na(markers_trim) & markers_trim != "NA"
                ]
            })
        })), # nolint
        MoreArgs = NULL
    ), .SDcols = c("cellMarker", "geneSymbol")]
}

cellmarker_download <- function(species) {
    if (is.null(cellmarker_database_external[[species]])) {
        data_link <- switch(species,
            human = "http://xteam.xbio.top/CellMarker/download/Human_cell_markers.txt",
            mouse = "http://xteam.xbio.top/CellMarker/download/Mouse_cell_markers.txt"
        )
        cli::cli_alert_info("Reading data from {.url {data_link}}")
        cellmarker_database_external[[species]] <- cellmarker_prepare(
            data.table::fread(data_link)
        )
        envir <- topenv(environment(NULL))
        unlockBinding("cellmarker_database_external", envir)
        assignInMyNamespace(
            "cellmarker_database_external",
            cellmarker_database_external
        )
        lockBinding("cellmarker_database_external", envir)
    }
    cellmarker_database_external[[species]]
}

cellmarker_database_external <- list(
    human = NULL, mouse = NULL
)
