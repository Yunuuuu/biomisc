#' Seach CellMarker database
#'
#' Search http://xteam.xbio.top/CellMarker
#' @param markers an atomic character, the markers to search in the CellMarker
#' database, can be the Gene Symbol, Gene ID, Protein Symbol or Protein ID
#' (usually starts with "P", "Q" or "O").
#' @param species a scalar string, "human" or "mouse".
#' @param internal logical value, indicates whether to use internal CellMarker
#' data. If `NULL`, this will be determined automatically; if the CellMarker
#' data has been downloaded once, namely, we have already used this function
#' once with a `internal` value `FALSE`, then the `NULL` will indicate `FALSE`.
#' Otherwise `TRUE`. The internal data was downloaded from CellMarker on
#' 2022-12-04.
#' @return a data.frame of the searching results. A column named `targeted`
#' containing the matched markers from CellMarker data, the row will be sorted
#' descendingly by the number of matched markers
#' @export
cellmarker_search <- function(markers, species = "human", internal = NULL) {
    # nolint start
    data <- data.table::copy(cellmarker_get(species, internal))
    data[, targeted := lapply(gene_list, function(.genes, .markers) {
        .genes[tolower(.genes) %in% .markers]
    }, .markers = tolower(markers))]
    data[, targeted_size := lengths(targeted)]
    data[, targeted_prop := targeted_size / length(markers)]
    data.table::setcolorder(
        data, c("targeted", "targeted_size", "targeted_prop"),
        after = "CellOntologyID"
    )
    data.table::setcolorder(
        data, intersect(cellmarker_gene_cols, names(data)),
        after = "gene_list"
    )
    data <- data[vapply(targeted, function(x) length(x) > 0L, logical(1L))][
        order(-targeted_size, na.last = TRUE)
    ]
    # nolint end
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
        function(...) {
            Reduce(union, list(...))
        },
        unname(lapply(.SD, function(markers) {
            markers_list <- strsplit(
                gsub("\\s*\\[\\s*|\\s*\\]\\s*", "", markers, perl = TRUE),
                "\\s*,\\s*",
                perl = TRUE
            )
            lapply(markers_list, function(markers_trim) {
                markers_trim <- trimws(markers_trim, "both", "[\\h\\v]")
                markers_trim[!is.na(markers_trim) & markers_trim != "NA"]
            })
        })),
        MoreArgs = NULL
    ), .SDcols = intersect(cellmarker_gene_cols, names(data))]
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
        utils::assignInMyNamespace(
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
cellmarker_gene_cols <- c(
    "cellMarker", "geneSymbol", "geneID", "proteinName", "proteinID"
)
