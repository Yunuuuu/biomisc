#' Identify consensus modules.
#' @param basis_list A list of sub-list (length >= 2) of basis matrix for each
#' NMF run.
#' @param method Specifying the method to define consensus modules. Defatils see
#' references.
#' @param min_size Minimal number of genes in a module.
#' @param ... Other argument passed to specific methods.
#'
#'  If method is "barkley", these can be follows:
#'   * min_jaccard: Minimal jaccard index to consider significant overlapped
#'     between two modules. Default: `0.05`.
#'   * s_min_jaccard: Modules significantly overlapped fewer than
#'     `s_min_jaccard` modules are removed. Default: `3L`.
#'   * cluster_fn: Function used to find communities (ususally function starts
#'     with "cluster_" in igraph package). Always accept a "graph" object.
#'     Default: `igraph::cluster_infomap`
#'   * v_min: Gene–gene connections were filtered out if they occurred in fewer
#'     than `v_min` individual tumor modules. Default: `2L`.
#'   * s_min: genes connected with fewer than `s_min` genes were removed.
#'     Default: `3L`.
#'
#'  If method is "gavish", these can be follows:
#'   * module_size: Expected size of genes in a module. Default: `50L`,
#'   * coef_threthold: Minimal coefficient  (defined by `basis`) of genes in a
#'     module.  Default: `0`.
#'   * min_intra_sim: Minimal intra-tumor overlap index, robust within the
#'     tumour (a program that is represented by several similar NMF programs, as
#'     defined for the same tumour when analysed by multiple NMF rank values;
#'     two NMF programs were considered as similar if they had at least 75% gene
#'     overlap defined by overlap). Default: `0.7`.
#'   * min_intra_size: Minimal number of significant programs (with jaccard
#'     index > min_intra_sim) to be considered robust within the tumor. If
#'     `NULL`, the number of NMF runs is used. Default: `NULL`.
#'   * min_inter_sim: Minimal inter-tumor jaccard index. robust across tumours
#'     (NMF programs that had at least 20% similarity (by top 50 genes) with any
#'     NMF program in any of the other tumours analysed). Default: `0.2`.
#'   * min_inter_size: Minimal number of significant programs (with overlap
#'     index > min_inter_sim) to be considered robust across tuomor. Default:
#'     `1`.
#'   * redundant_sim: Non-redundant within the tumour (within each tumour, NMF
#'     programs were ranked by their similarity (gene overlap) with NMFs from
#'     other tumours and selected in decreasing order; once an NMF was selected,
#'     any other NMF within the tumour that had 20% overlap (or more) with the
#'     selected NMF was removed, to avoid redundancy. Default: `0.2`.
#'   * min_consensus_sim: Minimal similarity index of gene intersection with
#'     other programs (across all tumor) to be considered significant. Default:
#'     `0.2`.
#'   * founder_intersection: Minimal size or fraction (if ranges from `0` to
#'     `1`) of significant intersected programs to define the first NMF program
#'     in a cluster. Default: `0.2`.
#'   * min_overlap_index: Minimal overlap cutoff for adding a new NMF program
#'     to the forming cluster. Default: `0.2`.
#'   * overlap_fn: Function to define overlap index between two gene sets.
#'     Default: `jaccard_index`.
#' @return A Consensus_module list.
#' @references
#' - Gavish, A., Tyler, M., Greenwald, A.C. et al. Hallmarks of transcriptional
#' intratumour heterogeneity across a thousand tumours. Nature 618, 598–606
#' (2023). https://doi.org/10.1038/s41586-023-06130-4
#' - Barkley, D., Moncada, R., Pour, M. et al. Cancer cell states recur across
#' tumor types and form specific interactions with the tumor microenvironment.
#' Nat Genet 54, 1192–1201 (2022). https://doi.org/10.1038/s41588-022-01141-9
#' @seealso
#' - <https://github.com/tiroshlab/3ca>
#' - <https://github.com/yanailab/PanCancer>
#' @export
nmf_consensus <- function(basis_list, method = c("barkley", "gavish"), min_size = 5L, ...) {
    method <- match.arg(method, c("barkley", "gavish"))
    assert_(
        basis_list, function(list) {
            is.list(list) && all(vapply(list, function(sublist) {
                is.list(sublist) && length(sublist) >= 2L &&
                    all(vapply(sublist, is.matrix, logical(1L)))
            }, logical(1L)))
        }, "a list of sublist (length >= 2) of {.cls matrix}"
    )
    # failed_elements <- NULL
    # for (i in seq_along(basis_list)) {
    #     if (!rlang::is_named(basis_list[[i]]) ||
    #         anyDuplicated(names(basis_list[[i]]))) {
    #         failed_elements <- c(failed_elements, i)
    #     }
    # }
    # if (length(failed_elements)) {
    #     if (!is.null(names(basis_list))) {
    #         failed_elements <- names(basis_list)[failed_elements]
    #     }
    #     cli::cli_abort(c(
    #         "All elements in {.arg module_list} must be named and shouldn't contain blank ({.val } and {.val NA_character_}) or duplicated names.", # nolint
    #         x = "Failed elements: {.val {failed_elements}}"
    #     )) # nolint
    # }
    switch(method,
        barkley = nmf_consensus_barkley(
            basis_list = basis_list,
            ..., min_size = min_size
        ),
        gavish = nmf_consensus_gavish(
            basis_list = basis_list,
            ..., min_size = min_size
        )
    )
}

#' @export
print.consensus_module <- function(x, ...) {
    for (at in setdiff(names(attributes(x)), "names")) {
        attr(x, at) <- NULL
    }
    print(x)
}

nmf_consensus_barkley <- function(basis_list, min_jaccard = 0.05, s_min_jaccard = 3L, v_min = 2L, s_min = 3L, cluster_fn = igraph::cluster_infomap, ..., min_size = 5L) {
    assert_pkg("igraph")
    module_list <- lapply(basis_list, nmf_modules_barkley, min_size = min_size)
    if (length(module_list) == 0L) {
        return(list())
    }
    assert_(cluster_fn, is.function, "a function")

    # https://github.com/yanailab/PanCancer/blob/49e7b270ec55dbe72076a5cae516ff0931fe7fe4/Finding.R#L349
    # retain modules with at least 5% overlap (by Jaccard index) with at least
    # two other modules
    raw_modules <- module_list
    names(module_list) <- seq_along(module_list)
    all_modules <- unlist(module_list, recursive = FALSE, use.names = TRUE)
    sim <- sapply(all_modules, function(x) {
        vapply(all_modules, jaccard_index, numeric(1L), x = x)
    })
    kept_modules <- rownames(sim)[rowSums(sim > min_jaccard) >= s_min_jaccard]
    if (length(kept_modules) == 0L) {
        cli::cli_abort("No modules to proceed")
    }
    all_modules <- all_modules[kept_modules]
    module_list <- mapply(
        function(x, i) {
            x[paste(i, names(x), sep = ".") %in% kept_modules]
        },
        x = module_list, i = names(module_list),
        SIMPLIFY = FALSE, USE.NAMES = TRUE
    )

    # Gene–gene connections ---------------------------------------
    # filtered out if they occurred in fewer than two individual tumor modules
    # ta <- table(unlist(all_modules, use.names = FALSE))
    # features <- names(ta)[ta > v_min]

    features <- unique(unlist(all_modules, use.names = FALSE))
    l <- length(features)
    adj_zero <- matrix(0L, nrow = l, ncol = l)
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
    modules <- igraph::communities(cluster_fn(g, ...))
    modules <- modules[lengths(modules) >= min_size]
    names(modules) <- paste0("MP", seq_along(modules))

    structure(
        modules,
        raw_programs = raw_modules,
        full_adj = full_adj,
        adj_list = adj_list,
        class = c("list", "consensus_module")
    )
}

nmf_modules_barkley <- function(basis_list, min_size = 5L) {
    # assert_(nmf_list, function(x) {
    #     is.list(x) && all(vapply(x, function(i) {
    #         inherits(i, "NMFfit")
    #     }, logical(1L)))
    # }, "a list of {.cls NMFfit}", cross_msg = NULL)

    # https://github.com/yanailab/PanCancer/blob/49e7b270ec55dbe72076a5cae516ff0931fe7fe4/part1.R#L67
    # ranks <- vapply(nmf_list, NMF::nbasis, integer(1L))
    ranks <- vapply(basis_list, ncol, integer(1L))

    # put NMFfit with the largest rank firstly
    # so we can use which.min to define the first near one with the largest rank
    rank_order <- order(ranks, decreasing = TRUE)
    ranks <- ranks[rank_order]
    basis_list <- basis_list[rank_order]
    modules_list <- lapply(basis_list, function(scores) {
        # Remove if fewer than min_size genes
        modules <- nmf_modules_internal(scores)
        scores <- scores[, lengths(modules) >= min_size, drop = FALSE]
        if (length(scores) == 0L) {
            return(character())
        }

        # Find modules
        modules <- nmf_modules_internal(scores)
        # names(modules) <- vapply(modules, `[[`, character(1L), 1L)
        names(modules) <- paste0("module", seq_along(modules))
        # names(modules) <- make.names(names(modules))
        modules
    })
    best_rank_idx <- which.min(ranks - lengths(modules_list))
    if (length(best_rank_idx)) {
        # nmf_out <- basis_list[[best_rank_idx]]
        modules <- modules_list[[best_rank_idx]]
    } else {
        cli::cli_abort("Cannot determine the best rank")
    }
    modules
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

jaccard_index <- function(x, y) {
    length(intersect(x, y)) / length(union(x, y))
}

nmf_consensus_gavish <- function(basis_list, min_size = 3L, coef_threthold = 0, module_size = 50L, min_intra_sim = 0.7, min_intra_size = NULL, min_inter_sim = 0.2, min_inter_size = 1L, redundant_sim = 0.2, min_consensus_sim = 0.2, founder_intersection = 0.2, min_overlap_index = 0.2, overlap_fn = jaccard_index) {
    assert_(overlap_fn, is.function, "a function")
    # Modified from <https://github.com/tiroshlab/3ca/blob/main/ITH_hallmarks/Generating_MPs/Generate_Meta_Programs.R>
    ### Parameters for clustering
    # min_intersection: the minimal intersection cutoff for defining the
    # first NMF program in a cluster;

    # min_overlap_index: the minimal intersection cutoff for adding a new
    # NMF to the forming cluster;

    # founder_intersection: the minimal group size to consider for defining
    # the first NMF program in a cluster

    # nmf programs -------------------------------------------
    nmf_program_list <- lapply(basis_list, function(rank_basis) {
        # use the integer index as the name
        names(rank_basis) <- seq_along(rank_basis)
        intra_program_list <- lapply(rank_basis, function(basis) {
            programs <- apply(basis, 2, function(y) {
                y <- y[y > coef_threthold]
                names(sort(y, decreasing = TRUE))[seq_len(min(length(y), module_size))]
            }, simplify = FALSE)
            # use the integer index as the name
            names(programs) <- seq_len(ncol(basis))
            programs[lengths(programs) >= min_size]
        })
        intra_program_list <- unlist(intra_program_list,
            recursive = FALSE, use.names = TRUE
        )
        l <- length(intra_program_list)
        if (l <= 1L) {
            list()
        } else {
            intra_similarity <- sapply(intra_program_list, function(x) {
                vapply(intra_program_list, function(y) {
                    overlap_fn(x, y)
                }, numeric(1L))
            })
            intra_program_list[
                rowSums(intra_similarity >= min_intra_sim) - 1L >=
                    (min_intra_size %||% (length(rank_basis) - 1L))
            ]
        }
    })
    # use the integer index as the name
    names(nmf_program_list) <- seq_along(nmf_program_list)

    recur_program_list <- lapply(names(nmf_program_list), function(id) {
        curr_program <- nmf_program_list[[id]]
        if (length(curr_program) == 0L) {
            return(list())
        }
        others <- nmf_program_list[setdiff(names(nmf_program_list), id)]
        others <- others[lengths(others) > 0L]
        inter_counts <- vapply(curr_program, function(program) {
            sum(vapply(others, function(other) {
                any(vapply(other, function(y) {
                    overlap_fn(program, y)
                }, numeric(1L)) >= min_inter_sim)
            }, numeric(1L)))
        }, numeric(1L))
        recur_idx <- inter_counts >= min_inter_size
        curr_program <- curr_program[recur_idx]
        inter_counts <- inter_counts[recur_idx]
        curr_program <- curr_program[order(inter_counts, decreasing = TRUE)]
        # filter redundant programs
        idx <- NULL
        for (i in seq_along(curr_program)) {
            if (length(idx) == 0L) {
                idx <- c(idx, i)
            } else {
                intra_similarity <- vapply(idx, function(i2) {
                    overlap_fn(curr_program[[i]], curr_program[[i2]])
                }, numeric(1L))
                if (max(intra_similarity) >= redundant_sim) {
                    next
                } else {
                    idx <- c(idx, i)
                }
            }
        }
        curr_program[idx]
    })
    names(recur_program_list) <- names(nmf_program_list)
    # names pattern: study_idx:nmf_rank_idx:nmf_factor_idx
    nmf_programs <- unlist(recur_program_list,
        recursive = FALSE, use.names = TRUE
    )

    # define consensus modules ---------------------------------
    # calculate similarity between programs
    sim_matrix <- sapply(nmf_programs, function(x) {
        vapply(nmf_programs, function(y) overlap_fn(x, y), numeric(1L))
    })
    sim_matrix_orig <- sim_matrix
    sorted_intersection <- sort(
        rowSums(sim_matrix >= min_consensus_sim) - 1L,
        decreasing = TRUE
    )
    program_list <- list() # Every entry contains the NMFs of a chosen cluster
    MP_list <- list()
    cur_program <- NULL
    k <- 1L

    if (founder_intersection > 0L && founder_intersection < 1L) {
        founder_intersection_threshold <- (ncol(sim_matrix) - 1L) *
            founder_intersection
    } else {
        founder_intersection_threshold <- founder_intersection
    }
    while (sorted_intersection[1] >= founder_intersection_threshold) {
        cur_program_id <- names(sorted_intersection[1])
        cur_program <- c(cur_program, cur_program_id)

        ### intersection between all remaining NMFs and Genes in MP
        genes_mp <- nmf_programs[[cur_program_id]]
        # Genes in the forming MP are first chosen to be those in the first NMF.
        # genes_mp always has only 50 genes and evolves during the formation of
        # the cluster
        nmf_programs <- nmf_programs[
            setdiff(names(nmf_programs), cur_program_id)
        ] # remove selected NMF
        if (length(nmf_programs) == 0L) {
            program_list[[paste0("MP", k)]] <- cur_program
            MP_list[[paste0("MP", k)]] <- genes_mp
            break
        }
        ovaerlap_genes_mp <- sapply(nmf_programs, function(y) {
            overlap_fn(genes_mp, y)
        })
        # intersection between all other NMFs and genes_mp
        ovaerlap_genes_mp <- sort(ovaerlap_genes_mp, decreasing = TRUE)
        all_genes_mp <- genes_mp
        # has genes in all NMFs in the current cluster, for redefining genes_mp
        # after adding a new NMF
        ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.

        while (ovaerlap_genes_mp[1] >= min_overlap_index) {
            next_program_id <- names(ovaerlap_genes_mp[1])
            cur_program <- c(cur_program, next_program_id)
            next_gene_mp <- nmf_programs[[next_program_id]]
            all_genes_mp <- c(all_genes_mp, next_gene_mp)
            all_genes_mp_counts <- sort(table(all_genes_mp), decreasing = TRUE)
            # genes_mp is newly defined each time according to all NMFs in the
            # current cluster
            if (length(all_genes_mp_counts) > module_size) {
                genes_at_border <- all_genes_mp_counts[
                    all_genes_mp_counts == all_genes_mp_counts[[module_size]]
                ] ### genes with overlap equal to the 50th gene
            } else {
                genes_at_border <- NULL
            }

            if (length(genes_at_border) > 1L) {
                ### Sort last genes in genes_at_border according to maximal NMF gene scores
                ### Run across all NMF programs in cur_program and extract NMF scores for each gene
                nmf_coef <- NULL
                for (i in cur_program) {
                    curr_study <- strsplit(i, ".", fixed = TRUE)[[1]]
                    curr_study <- as.integer(curr_study)
                    curr_basis <- basis_list[[curr_study[1:2]]]
                    curr_nmf_coef <- curr_basis[
                        intersect(names(genes_at_border), rownames(curr_basis)),
                        curr_study[[3L]],
                        drop = TRUE
                    ]
                    nmf_coef <- c(nmf_coef, curr_nmf_coef)
                }
                nmf_coef <- sort(nmf_coef, decreasing = TRUE)
                nmf_coef <- nmf_coef[unique(names(nmf_coef))]
                genes_mp <- c(
                    names(all_genes_mp_counts[
                        all_genes_mp_counts > all_genes_mp_counts[
                            min(module_size, length(all_genes_mp_counts))
                        ]
                    ]),
                    names(nmf_coef)
                )
            } else {
                genes_mp <- names(all_genes_mp_counts)
            }
            genes_mp <- genes_mp[seq_len(min(module_size, length(genes_mp)))]
            remain_programs <- setdiff(names(nmf_programs), next_program_id)
            nmf_programs <- nmf_programs[remain_programs] # remove selected NMF
            if (length(nmf_programs) == 0L) break
            ovaerlap_genes_mp <- sapply(nmf_programs, function(y) {
                overlap_fn(genes_mp, y)
            })
            # intersection between all other NMFs and genes_mp
            ovaerlap_genes_mp <- sort(ovaerlap_genes_mp, decreasing = TRUE)
        }
        program_list[[paste0("MP", k)]] <- cur_program
        MP_list[[paste0("MP", k)]] <- genes_mp
        sim_matrix <- sim_matrix[
            setdiff(rownames(sim_matrix), cur_program),
            setdiff(colnames(sim_matrix), cur_program),
            drop = FALSE
        ] # Remove current chosen cluster
        if (length(sim_matrix) == 0L) break
        sorted_intersection <- sort(
            rowSums(sim_matrix >= min_consensus_sim) - 1L,
            decreasing = TRUE
        )
        if (founder_intersection > 0L && founder_intersection < 1L) {
            founder_intersection_threshold <- (ncol(sim_matrix) - 1L) *
                founder_intersection
        } else {
            founder_intersection_threshold <- founder_intersection
        }
        cur_program <- NULL
        k <- k + 1
    }

    # order final program data based on clustering output --------
    # hierarchical clustering of the similarity matrix
    nmf_hc <- stats::hclust(stats::as.dist(1 - sim_matrix_orig),
        method = "average"
    )
    nmf_hc <- stats::reorder(
        stats::as.dendrogram(nmf_hc),
        colMeans(sim_matrix_orig)
    )
    order_hc <- stats::order.dendrogram(nmf_hc)
    sim_matrix_orig <- sim_matrix_orig[order_hc, order_hc]
    idxs_sorted <- unlist(program_list, recursive = FALSE, use.names = FALSE)
    sim_matrix_group <- rep(names(program_list), times = lengths(program_list))
    non_classified <- setdiff(colnames(sim_matrix_orig), idxs_sorted)
    idxs_sorted <- c(idxs_sorted, non_classified)
    sim_matrix_group <- c(
        sim_matrix_group,
        rep_len("none", length(non_classified))
    )
    sim_matrix_group <- factor(sim_matrix_group, unique(sim_matrix_group))
    sim_matrix_orig <- sim_matrix_orig[idxs_sorted, idxs_sorted]

    recur_program_list <- mapply(
        function(programs, basis) {
            if (length(programs) == 0L) {
                return(list())
            }
            names(programs) <- nmf_consensus_parse_name(
                names(programs),
                basis = basis
            )
            programs
        },
        basis = basis_list,
        programs = recur_program_list,
        SIMPLIFY = FALSE, USE.NAMES = TRUE
    )
    program_list <- lapply(program_list, function(programs) {
        nmf_consensus_parse_name(programs, basis_list)
    })
    rownames(sim_matrix_orig) <- nmf_consensus_parse_name(
        rownames(sim_matrix_orig), basis_list
    )
    colnames(sim_matrix_orig) <- nmf_consensus_parse_name(
        colnames(sim_matrix_orig), basis_list
    )
    structure(
        MP_list,
        elements = program_list,
        raw_programs = recur_program_list,
        sim_matrix = sim_matrix_orig,
        sim_matrix_group = sim_matrix_group,
        class = c("list", "consensus_module")
    )
}

nmf_consensus_parse_name <- function(idx_string, basis_list = NULL, basis = NULL) {
    vapply(strsplit(idx_string, ".", fixed = TRUE), function(x) {
        if (is.null(basis_list)) {
            all_names <- nmf_consensus_parse_name_rank(as.integer(x), basis)
        } else {
            all_names <- nmf_consensus_parse_name_project(
                as.integer(x), basis_list
            )
        }
        paste(all_names, collapse = ".")
    }, character(1L))
}

nmf_consensus_parse_name_project <- function(idx, basis_list) {
    project_id <- (names(basis_list) %||% seq_along(basis_list))[
        idx[[1L]]
    ]
    basis <- basis_list[[idx[[1L]]]]
    c(project_id, nmf_consensus_parse_name_rank(idx[-1L], basis))
}

nmf_consensus_parse_name_rank <- function(idx, basis) {
    rank_id <- (names(basis) %||% seq_along(basis))[idx[1L]]
    curr_basis <- basis[[idx[[1L]]]]
    c(rank_id, nmf_consensus_parse_name_factor(idx[[2L]], curr_basis))
}

nmf_consensus_parse_name_factor <- function(idx, basis) {
    (colnames(basis) %||% seq_len(ncol(basis)))[idx]
}
