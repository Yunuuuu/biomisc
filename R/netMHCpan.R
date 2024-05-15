#' Parse NetMHCpan result
#' @param file A string of the result from NetMHCpan.
#' @param ids A character of identifier of the fasta sequence which was input
#' into NetMHCpan.
#' @param rth Rank Threshold for high binding peptides.
#' @param rlt Rank Threshold for low binding peptides.
#' @details
#'  BindLevel: (SB: strong binder, WB: weak binder). The peptide will be
#'  identified as a strong binder if the % Rank is below the specified threshold
#'  for the strong binders (`rth`). The peptide will be identified as a weak
#'  binder if the % Rank is above the threshold of the strong binders (`rth`)
#'  but below the specified threshold for the weak binders (`rlt`).
#' @return A [data.table][data.table::data.table].
#' @seealso <https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.1/>
#' @export
parse_netmhcpan <- function(file, ids = NULL, rth = 0.5, rlt = 2.0) {
    data <- data.table::fread(file, header = FALSE, skip = 2L)

    nms <- unlist(
        data.table::fread(file, header = FALSE, nrows = 1L, skip = 1L),
        recursive = FALSE, use.names = FALSE
    )
    hla_genotype <- unlist(data.table::fread(file, header = FALSE, nrows = 1L),
        recursive = FALSE, use.names = FALSE
    )
    hla_genotype <- hla_genotype[!is.na(hla_genotype)]
    measure_vars <- paste(
        rep(hla_genotype, each = 6L),
        nms[4:(length(nms) - 2L)],
        sep = "::"
    )
    nms[4:(length(nms) - 2L)] <- measure_vars
    # nolint start
    data[, identifiers := cumsum(
        c(0L, diff(nchar(Peptide))) > 0L |
            V3 != data.table::shift(V3, type = "lag", fill = head(V3, n = 1L))
    ) + 1L]
    data.table::setnames(data, nms)
    if (!is.null(ids)) {
        data[, identifiers := ids[identifiers]]
    }
    data <- data.table::melt(
        data,
        id.vars = c(names(data)[1:3], "identifiers"),
        measure.vars = split(
            measure_vars,
            sub("^.+::(.+)$", "\\1", measure_vars)
        ),
        variable.factor = FALSE, variable.name = "HLA"
    )
    data[, HLA := hla_genotype[as.integer(HLA)]]
    data[, BindLevel := data.table::fcase(
        EL_Rank < rth, "SB",
        EL_Rank < rlt, "WB"
    )][]
    # nolint end
}

utils::globalVariables(
    c("identifiers", "Peptide", "HLA", "BindLevel", "EL_Rank")
)
