#' Identify modules across different NMF ranks.
#'
#' This function aims to ascertain the optimal rank for a given NMF object list
#' and subsequently identify the modules linked to the aforementioned
#' determination.
#' @param nmf_list A named list of NMF (identified by [NMF][NMF::NMF]).
#' @param threshold Indicates the minimum number of genes required to filter a
#'  module, failing which, it shall not be considered.
#' @return A list, including a curated selection of NMF objects and the idenfied
#'  modules.
#' @references
#' Barkley, D., Moncada, R., Pour, M. et al. Cancer cell states recur across
#' tumor types and form specific interactions with the tumor microenvironment.
#' Nat Genet 54, 1192–1201 (2022). https://doi.org/10.1038/s41588-022-01141-9
#' @seealso
#' <https://github.com/yanailab/PanCancer>
#' @export
nmf_modules <- function(nmf_list, threshold = 5L) {
    assert_pkg("NMF")
    assert_class(nmf_list, function(x) {
        is.list(x) && all(vapply(x, function(i) {
            inherits(i, "NMFfit")
        }, logical(1L)))
    }, "{.cls list} of {.cls NMFfit}", cross_msg = NULL)

    # https://github.com/yanailab/PanCancer/blob/49e7b270ec55dbe72076a5cae516ff0931fe7fe4/part1.R#L67
    ranks <- vapply(nmf_list, NMF::nbasis, integer(1L))

    # put NMFfit with the largest rank firstly
    # so we can use which.min to define the first near one with the largest rank
    rank_order <- order(ranks, decreasing = TRUE)
    ranks <- ranks[rank_order]
    nmf_list <- nmf_list[rank_order]
    modules_list <- lapply(nmf_list, function(nmf_out) {
        scores <- NMF::basis(nmf_out)

        # Remove if fewer than threshold genes
        modules <- nmf_modules_internal(scores)
        scores <- scores[, lengths(modules) >= threshold, drop = FALSE]
        if (length(scores) == 0L) {
            return(character())
        }

        # Find modules
        modules <- nmf_modules_internal(scores)
        # names(modules) <- vapply(modules, `[[`, character(1L), 1L)
        names(modules) <- paste("module", seq_along(modules), sep = "_")
        # names(modules) <- make.names(names(modules))
        modules
    })
    best_rank_idx <- which.min(ranks - lengths(modules_list))
    if (length(best_rank_idx)) {
        nmf_out <- nmf_list[[best_rank_idx]]
        modules <- modules_list[[best_rank_idx]]
        class(modules) <- "NMFmodule"
    } else {
        cli::cli_abort("Cannot determine the best rank")
    }
    list(nmf = nmf_out, modules = modules)
}

#' @export
print.NMFmodule <- function(x, ...) {
    print(unclass(x))
}

nmf_modules_internal <- function(scores) {
    contributions <- -t(t(scores) / colMeans(scores))
    ranks_x <- t(apply(contributions, 1L, rank))
    ranks_y <- apply(contributions, 2L, rank)
    # test_ranks_y <- ranks_y
    # for (i in seq_len(ncol(scores))) {
    #     ranks_y[ranks_x[, i] > 1, i] <- Inf
    # }
    # identical(ranks_y, test_ranks_y)
    ranks_y[ranks_x > 1L] <- Inf
    apply(ranks_y, 2L, function(m) {
        a <- sort(m[is.finite(m)])
        names(a[a == seq_along(a)])
    }, simplify = FALSE)
}
