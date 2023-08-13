#' Identify consensus modules.
#' @param module_list A named list (defining tumor type) of named sub-list of
#'  atomic characters (each stands for the identified modules in specific
#'  tumor). The names of the list must be unique and blank names are not
#'  allowed.
#' @param cluster_fn Function used to find communities (ususally function starts
#'  with "cluster_" in igraph package). Always accept a "graph" object.
#' @param ... Other argument passed to `cluster_fn`.
#' @param v_min Gene–gene connections were filtered out if they occurred in
#' fewer than `v_min` individual tumor modules.
#' @param s_min genes connected with fewer than `s_min` genes were removed.
#' @param threshold Indicates the minimum number of genes required to filter a
#'  module, failing which, it shall not be considered.
#' @return A Consensus_module list of filtered communities by `threshold`, each
#' identified by their members, with `raw_modules`, `full_adj` and `adj_list`
#' attached in the attributes.
#' @references
#' Barkley, D., Moncada, R., Pour, M. et al. Cancer cell states recur across
#' tumor types and form specific interactions with the tumor microenvironment.
#' Nat Genet 54, 1192–1201 (2022). https://doi.org/10.1038/s41588-022-01141-9
#' @seealso
#' <https://github.com/yanailab/PanCancer>
#' @export
consensus_modules <- function(module_list, cluster_fn = igraph::cluster_infomap, ..., v_min = 2L, s_min = 3L, threshold = 5L) {
    assert_pkg("igraph")
    assert_class(module_list, function(x) {
        is.list(x) && all(vapply(x, function(i) {
            is.list(i) && all(vapply(i, is.character, logical(1L)))
        }, logical(1L)))
    }, "{.cls list} of sub-list of atomic {.cls character}", cross_msg = NULL)
    if (length(module_list) == 0L) {
        return(list())
    }
    if (!rlang::is_named(module_list) || anyDuplicated(names(module_list))) {
        cli::cli_abort("{.arg module_list} must be named and shouldn't contain blank ({.val } and {.val NA_character_}) or duplicated names.") # nolint
    }
    failed_elements <- NULL
    for (i in names(module_list)) {
        if (!rlang::is_named(module_list[[i]]) ||
            anyDuplicated(names(module_list[[i]]))) {
            failed_elements <- c(failed_elements, i)
        }
    }
    if (length(failed_elements)) {
        cli::cli_abort(c(
            "All elements in {.arg module_list} must be named and shouldn't contain blank ({.val } and {.val NA_character_}) or duplicated names.", # nolint
            x = "Failed elements: {.val {failed_elements}}"
        )) # nolint
    }
    assert_class(cluster_fn, is.function, "{.cls function}", cross_msg = NULL)
    # https://github.com/yanailab/PanCancer/blob/49e7b270ec55dbe72076a5cae516ff0931fe7fe4/Finding.R#L349
    # retain modules with at least 5% overlap (by Jaccard index) with at least
    # two other modules
    all <- unlist(module_list, recursive = FALSE, use.names = TRUE)
    sim <- sapply(all, function(x) {
        vapply(all, jaccard_index, numeric(1L), x = x)
    })
    kept_modules <- rownames(sim)[rowSums(sim > 0.05) >= 3L]
    if (length(kept_modules) == 0L) {
        cli::cli_abort("No modules to proceed")
    }
    all <- all[kept_modules]
    module_list <- mapply(
        function(x, i) {
            x[paste(i, names(x), sep = ".") %in% kept_modules]
        },
        x = module_list, i = names(module_list),
        SIMPLIFY = FALSE, USE.NAMES = TRUE
    )

    # Gene–gene connections ---------------------------------------
    # filtered out if they occurred in fewer than two individual tumor modules
    # ta <- table(unlist(all, use.names = FALSE))
    # features <- names(ta)[ta > v_min]

    features <- unique(unlist(all, use.names = FALSE))
    adj_zero <- matrix(0L, nrow = length(features), ncol = length(features))
    rownames(adj_zero) <- colnames(adj_zero) <- features

    # Adjacency matrix, list by cancer
    adj_list <- lapply(module_list, function(modules, init = adj_zero) {
        for (i in seq_along(modules)) {
            mod <- intersect(modules[[i]], features)
            init[mod, mod] <- init[mod, mod] + 1L
        }
        diag(init) <- 0L
        init
    })
    full_adj <- Reduce(`+`, adj_list)

    # Remove low connections
    # Gene–gene connections were filtered out if they occurred in fewer than two
    # individual tumor modules, and genes with fewer than three connections were
    # removed.
    full_adj[full_adj < v_min] <- 0L
    is_keep <- rowSums(full_adj > 0L) >= s_min
    full_adj <- full_adj[is_keep, is_keep]

    # Cluster ---------------------------------------------------
    g <- igraph::graph_from_adjacency_matrix(full_adj,
        diag = FALSE, mode = "undirected", weighted = TRUE
    )
    raw_modules <- igraph::communities(cluster_fn(graph = g, ...))
    modules <- raw_modules[lengths(raw_modules) >= threshold]
    names(modules) <- paste0("consensus_module_", seq_along(modules))
    structure(modules,
        raw_modules = raw_modules,
        full_adj = full_adj,
        adj_list = adj_list,
        class = c("list", "consensus_module")
    )
}

#' @export
print.consensus_module <- function(x, ...) {
    x <- unclass(x)
    attr(x, "raw_modules") <- NULL
    attr(x, "adj_list") <- NULL
    attr(x, "full_adj") <- NULL
    print(x)
}

jaccard_index <- function(x, y) {
    length(intersect(x, y)) / length(union(x, y))
}
