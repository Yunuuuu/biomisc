#' Identify consensus modules.
#' @param basis_list A list of sub-list of matrix for each NMF run.
#' @param module_size Expected size of genes in a module.
#' @param min_contribution Minimal contricution percent (defined by `basis`) of
#' genes in a module.
#' @param min_size Minimal number of genes in a module.
#' @param min_intra_sim Minimal intra-tumor jaccard index, robust within the
#' tumour (a program that is represented by several similar NMF programs, as
#' defined for the same tumour when analysed by multiple NMF rank values; two
#' NMF programs were considered as similar if they had at least 75% gene overlap
#' defined by jaccard index)
#' @param min_intra_size Minimal number of significant programs (with jaccard
#' index > min_intra_sim) to be considered robust within the tumor. If NULL, the
#' number of NMF runs is used.
#' @param min_inter_sim Minimal inter-tumor jaccard index. robust across tumours
#' (NMF programs that had at least 20% similarity (by top 50 genes) with any NMF
#' program in any of the other tumours analysed).
#' @param min_inter_size Minimal number of significant programs (with jaccard
#' index > min_inter_sim) to be considered robust across tuomor.
#' @param redundant_sim Non-redundant within the tumour (within each tumour, NMF
#' programs were ranked by their similarity (gene overlap) with NMFs from other
#' tumours and selected in decreasing order; once an NMF was selected, any other
#' NMF within the tumour that had 20% overlap (or more) with the selected NMF
#' was removed, to avoid redundancy.
#' @param min_intersection Minimal size of gene intersection with other programs
#' (across all tumor) to be considered significant.
#' @param founder_intersection_size Minimal size of significant intersected
#' programs to define the first NMF program in a cluster.
#' @param child_intersection_size Minimal intersection cutoff for adding a new
#' NMF program to the forming cluster.
#' @references
#' Gavish, A., Tyler, M., Greenwald, A.C. et al. Hallmarks of transcriptional
#' intratumour heterogeneity across a thousand tumours. Nature 618, 598–606
#' (2023). https://doi.org/10.1038/s41586-023-06130-4
#' @seealso
#' <https://github.com/tiroshlab/3ca>
#' @export
nmf_consensus <- function(basis_list, module_size = 50L, min_contribution = 0.02, min_size = 3L, min_intra_sim = 0.7, min_intra_size = NULL, min_inter_sim = 0.2, min_inter_size = 2L, redundant_sim = 0.2, min_intersection = 0.2, founder_intersection_size = 0.2, child_intersection_size = 0.2) {
    assert_class(
        basis_list, function(list) {
            is.list(list) && all(vapply(list, function(sublist) {
                is.list(sublist) && length(sublist) >= 2L &&
                    all(vapply(sublist, is.matrix, logical(1L)))
            }, logical(1L)))
        }, "{.cls list} of sublist (length >= 2) of {.cls matrix}",
        cross_msg = NULL
    )
    # Modified from <https://github.com/tiroshlab/3ca/blob/main/ITH_hallmarks/Generating_MPs/Generate_Meta_Programs.R>
    ### Parameters for clustering
    # min_intersection: the minimal intersection cutoff for defining the
    # first NMF program in a cluster;

    # child_intersection_size: the minimal intersection cutoff for adding a new
    # NMF to the forming cluster;

    # founder_intersection_size: the minimal group size to consider for defining
    # the first NMF program in a cluster

    # nmf programs -------------------------------------------
    nmf_program_list <- lapply(basis_list, function(rank_basis) {
        # use the integer index as the name
        names(rank_basis) <- seq_along(rank_basis)
        intra_program_list <- lapply(rank_basis, function(basis) {
            programs <- apply(basis, 2, function(y) {
                sum_val <- sum(y)
                if (sum_val == 0L) {
                    return(character())
                }
                y <- y[y / sum_val >= min_contribution]
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
                sapply(intra_program_list, jaccard_index, x = x)
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
                    jaccard_index(program, y)
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
                    jaccard_index(curr_program[[i]], curr_program[[i2]])
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
    nmf_intersect <- sapply(nmf_programs, function(x) {
        sapply(nmf_programs, function(y) length(intersect(x, y)))
    })
    nmf_union <- sapply(nmf_programs, function(x) {
        sapply(nmf_programs, function(y) length(union(x, y)))
    })
    min_length <- sapply(nmf_programs, function(x) {
        sapply(nmf_programs, function(y) min(lengths(list(x, y))))
    })
    # hierarchical clustering of the similarity matrix
    nmf_intersect_hc <- stats::hclust(
        stats::as.dist(1 - (nmf_intersect / nmf_union)),
        method = "average"
    )
    nmf_intersect_hc <- stats::reorder(
        stats::as.dendrogram(nmf_intersect_hc),
        colMeans(nmf_intersect)
    )
    order_hc <- stats::order.dendrogram(nmf_intersect_hc)
    nmf_intersect <- nmf_intersect[order_hc, order_hc]
    min_length <- min_length[order_hc, order_hc]
    if (min_intersection > 0L && min_intersection < 1L) {
        min_intersection <- min_length * min_intersection
    }
    sorted_intersection <- sort(
        rowSums(nmf_intersect >= min_intersection) - 1L,
        decreasing = TRUE
    )
    program_list <- list() # Every entry contains the NMFs of a chosen cluster
    MP_list <- list()
    cur_program <- NULL
    k <- 1L

    if (founder_intersection_size > 0L && founder_intersection_size < 1L) {
        founder_intersection_threshold <- ncol(nmf_intersect) *
            founder_intersection_size
    } else {
        founder_intersection_threshold <- founder_intersection_size
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

        intersection_with_genes_mp <- nmf_intersect[
            names(nmf_programs), cur_program_id,
            drop = TRUE
        ]
        intersection_with_genes_mp <- sort(
            intersection_with_genes_mp,
            decreasing = TRUE
        ) # intersection between all other NMFs and genes_mp
        all_genes_mp <- genes_mp
        # has genes in all NMFs in the current cluster, for redefining genes_mp
        # after adding a new NMF
        ### Create gene list is composed of intersecting genes (in descending order by frequency). When the number of genes with a given frequency span bewond the 50th genes, they are sorted according to their NMF score.

        if (child_intersection_size > 0L && child_intersection_size < 1L) {
            child_intersection_threshold <- module_size * child_intersection_size
        } else {
            child_intersection_threshold <- child_intersection_size
        }
        child_intersection_threshold <- as.integer(child_intersection_threshold)
        while (intersection_with_genes_mp[1] >= child_intersection_threshold) {
            next_program_id <- names(intersection_with_genes_mp[1])
            cur_program <- c(cur_program, next_program_id)

            all_genes_mp <- c(all_genes_mp, nmf_programs[[next_program_id]])
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
                    module_basis <- basis_list[[curr_study[1:2]]]
                    curr_nmf_coef <- module_basis[
                        intersect(
                            names(genes_at_border), rownames(module_basis)
                        ),
                        curr_study[[3L]],
                        drop = TRUE
                    ]
                    nmf_coef <- c(nmf_coef, curr_nmf_coef)
                }
                nmf_coef <- sort(nmf_coef, decreasing = TRUE)
                nmf_coef <- nmf_coef[unique(names(nmf_coef))]
                all_genes_mp_counts <- c(
                    names(all_genes_mp_counts[
                        all_genes_mp_counts > all_genes_mp_counts[
                            min(module_size, length(all_genes_mp_counts))
                        ]
                    ]),
                    names(nmf_coef)
                )
            } else {
                all_genes_mp_counts <- names(all_genes_mp_counts)[
                    seq_len(min(module_size, length(all_genes_mp_counts)))
                ]
            }
            genes_mp <- all_genes_mp_counts[
                seq_len(min(module_size, length(all_genes_mp_counts)))
            ]
            remain_programs <- setdiff(names(nmf_programs), next_program_id)
            if (length(remain_programs) == 0L) break
            nmf_programs <- nmf_programs[remain_programs] # remove selected NMF
            intersection_with_genes_mp <- nmf_intersect[
                names(nmf_programs), next_program_id,
                drop = TRUE
            ]
            names(intersection_with_genes_mp) <- names(nmf_programs)
            intersection_with_genes_mp <- sort(
                intersection_with_genes_mp,
                decreasing = TRUE
            ) # intersection between all other NMFs and genes_mp
        }
        program_list[[paste0("MP", k)]] <- vapply(
            strsplit(cur_program, ".", fixed = TRUE), function(x) {
                idx <- as.integer(x)
                first_id <- (names(basis_list) %||% seq_along(basis_list))[
                    idx[[1L]]
                ]
                basis <- basis_list[[idx[[1L]]]]
                second_id <- (names(basis) %||% seq_along(basis))[idx[2L]]
                cur_basis <- basis[[idx[2]]]
                third_id <- (colnames(cur_basis) %||% seq_len(ncol(cur_basis)))[
                    idx[[3L]]
                ]
                paste(first_id, second_id, third_id, sep = "-")
            }, character(1L)
        )
        MP_list[[paste0("MP", k)]] <- genes_mp
        if (!is.null(min_intersection) && all(dim(min_intersection) == dim(nmf_intersect))) {
            min_intersection <- min_intersection[
                setdiff(rownames(min_intersection), cur_program),
                setdiff(colnames(min_intersection), cur_program),
                drop = FALSE
            ]
        }
        nmf_intersect <- nmf_intersect[
            setdiff(rownames(nmf_intersect), cur_program),
            setdiff(colnames(nmf_intersect), cur_program),
            drop = FALSE
        ] # Remove current chosen cluster
        if (length(nmf_intersect) == 0L) break
        sorted_intersection <- sort(
            rowSums(nmf_intersect >= min_intersection) - 1L,
            decreasing = TRUE
        )
        cur_program <- NULL
        k <- k + 1
    }
    recur_program_list <- mapply(
        function(programs, basis) {
            if (length(programs) == 0L) {
                return(list())
            }
            idx <- strsplit(names(programs), ".", fixed = TRUE)
            idx <- data.table::transpose(idx) # rank_id, factor_id
            rank_ids <- names(basis) %||% seq_along(basis)
            names(programs) <- mapply(function(x, y) {
                x <- as.integer(x)
                y <- as.integer(y)
                rank_id <- rank_ids[[x]]
                cur_basis <- basis[[x]]
                factor_ids <- colnames(cur_basis) %||% seq_len(ncol(cur_basis))
                factor_id <- factor_ids[[y]]
                paste(rank_id, factor_id, sep = "-")
            }, x = idx[[1L]], y = idx[[2L]], USE.NAMES = FALSE)
            programs
        },
        basis = basis_list,
        programs = recur_program_list,
        SIMPLIFY = FALSE, USE.NAMES = TRUE
    )
    list(
        recur_program_list = recur_program_list,
        meta_program = MP_list, clusters = program_list
    )
}
