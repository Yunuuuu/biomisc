## Rank Rank Hypergeometric Overlap
## based on Plaisier et al., Nucleic Acids Research, 2010
## list 1 is a data.frame from experiment 1 with two columns,
## column 1 is the Gene Identifier,
## column 2 is the signed ranking value (e.g. signed -log(p-value)
##        or fold change)
##
## list 2 is a data.frame from experiment 2 with two columns,
## column 1 is the Gene Identifier,
## column 2 is the signed ranking value (e.g. signed -log10(p-value)
##    or fold change)
## stepsize indicates how many genes to increase by
##    in each algorithm iteration
### Testing:
# n <- 112
# sample1 <- sample(n)
# sample2 <- sample(n)
# names(sample1) <- names(sample2) <- paste0("gene", seq_len(n))
# res1 <- calculate_hyper_overlap(
#     sample1, sample2, n, 10, "enrichment"
# )
# res2 <- RRHO:::numericListOverlap(
#     names(sample1), names(sample2), 10,
#     alternative = "enrichment"
# )
# all(dplyr::near(
#     res1$counts, res2$counts
# ))
# all(dplyr::near(
#     -log(res1$pvalue), res2$log.pval
# ))
# all(dplyr::near(
#     res1$counts, res2$counts
# ))
# res3 <- run_rrho(sample1, sample2, 1)
# set.seed(1)
# rrho_corrected_pval(res3, "perm", 10)
# rrho_corrected_pval(res3, "perm", 10)
# set.seed(1)
# rrho_corrected_pval(res3, "perm", 10)
# rrho_corrected_pval(res3, "perm", 10)
# res4 <- RRHO::RRHO(
#     as.data.frame(tibble::enframe(sample1)),
#     as.data.frame(tibble::enframe(sample2)),
#     1L,
#     alternative = "enrichment", plots = TRUE,
#     labels = c("sample1", "sample2"),
#     outputdir = tempdir()
# )
# all(dplyr::near(
#     res3$hypermat,
#     res4$hypermat
# ))
# all(dplyr::near(
#     res3$hypermat_counts,
#     res4$hypermat.counts
# ))
# all(dplyr::near(
#     res3$hypermat_signs,
#     res4$hypermat.signs
# ))
# upgenes1 <- rrho_sig_terms(res3)
# jet.colors <- colorRampPalette(
#     c(
#         "#00007F", "blue", "#007FFF", "cyan",
#         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"
#     )
# )
# # layout(matrix(c(rep(1, 5), 2), 1, 6, byrow = TRUE))
# image(res3$hypermat * res3$hypermat_signs,
#     xlab = "", ylab = "", col = jet.colors(100),
#     axes = FALSE, main = "Rank Rank Hypergeometric Overlap Map"
# )
# image(res4$hypermat * res4$hypermat.signs,
#     xlab = "", ylab = "", col = jet.colors(100),
#     axes = FALSE, main = "Rank Rank Hypergeometric Overlap Map"
# )
# rrho_heatmap(res3)
# rrho_corrected_pval(res3)
# rrho_corrected_pval(res3, "permutation", 20)
# temp_rrho_res <- readRDS("data-raw/diff_rrho_res.rds")
# rrho_heatmap(temp_rrho_res)
#' Rank-Rank Hypergeometric Overlap Test
#'
#' The function tests for significant overlap between two sorted lists using the method in the reference.
#'
#' @param list1 a named numeric vector, For differential gene expression, values
#' are often `-log10(P-value) * sign(effect)`.
#' @param list2 Same as list1
#' @param stepsize Controls the resolution of the test: how many items between
#' any two overlap tests.
#' @param alternative one of c("enrichment", "two.sided")
#' @param log_base Normally, `hypermat` in the results are
#' transformed by logarithm, this control the logarithm base. Just like the
#' `base` parameter in [base::log] function.
#' @export
run_rrho <- function(list1, list2, stepsize = NULL, alternative = NULL,
                     log_base = exp(1L)) {
    alternative <- match.arg(alternative, c("enrichment", "two.sided"))
    # remove NA value and zero value
    list1 <- list1[
        !is.na(list1) & !abs(list1 - 0L) < sqrt(.Machine$double.eps)
    ]
    list2 <- list2[
        !is.na(list2) & !abs(list2 - 0L) < sqrt(.Machine$double.eps)
    ]
    # keep items in both lists
    common_names <- intersect(
        names(list1), names(list2)
    )
    list1 <- list1[common_names]
    list2 <- list2[common_names]
    rlang::inform(
        sprintf("Found %d genes shared by both list.", length(common_names))
    )
    if (is.null(stepsize)) stepsize <- set_stepsize(list1, list2)

    ## Order both lists
    list1 <- sort(list1, decreasing = TRUE)
    list2 <- sort(list2, decreasing = TRUE)
    hyper_res <- calculate_hyper_overlap(
        list1, list2,
        n = length(common_names),
        stepsize = stepsize, alternative = alternative
    )

    hyperlogp <- -log(hyper_res$pvalue, base = log_base)
    list(
        hypermat = hyperlogp,
        hypermat_counts = hyper_res$counts,
        hypermat_signs = hyper_res$signs,
        hypermat_pvalue = hyper_res$pvalue,
        ranked_list = list(
            list1 = list1,
            list2 = list2
        ),
        items = common_names,
        nitems = length(common_names),
        stepsize = stepsize,
        log_base = log_base,
        alternative = alternative
    )
}


## Suggest default step size
set_stepsize <- function(list1, list2) {
    n1 <- length(list1)
    n2 <- length(list2)
    ceiling(sqrt(min(c(n1, n2))))
}

## Compute the overlaps between two *numeric* lists:
hyper_test <- function(sample1, sample2, alternative, n, tol) {
    count <- length(intersect(sample1, sample2))
    a <- length(sample1)
    b <- length(sample2)
    if (identical(alternative, "enrichment")) {
        pvalue <- stats::phyper(
            q = count - 1L, m = a, n = n - a + 1L,
            k = b, lower.tail = FALSE, log.p = FALSE
        )
        sign <- 1L
    } else if (identical(alternative, "two.sided")) {
        the_mean <- a * b / n
        sign <- as.integer(sign(count - the_mean))
        if (sign < 0L) {
            lower <- count
            upper <- 2L * the_mean - count
        } else {
            lower <- 2L * the_mean - count
            upper <- count
        }
        pval1 <- stats::phyper(
            q = lower + tol, m = a,
            n = n - a + 1L, k = b,
            lower.tail = TRUE
        )
        pval2 <- stats::phyper(
            q = upper - tol, m = a,
            n = n - a + 1L, k = b,
            lower.tail = FALSE
        )
        pvalue <- pval1 + pval2
    }
    c(
        count = count,
        pvalue = pvalue,
        sign = sign
    )
}

calculate_hyper_overlap <- function(list1, list2, n, stepsize, alternative, tol = 0.5) {
    col_ids <- row_ids <- seq(1L, n, by = stepsize)
    indexes <- expand.grid(
        row_ids = row_ids,
        col_ids = col_ids
    )
    indexes <- as.matrix(indexes)
    overlaps <- apply(indexes, 1L, function(x) {
        hyper_test(
            names(list1)[seq_len(x[["row_ids"]])],
            names(list2)[seq_len(x[["col_ids"]])],
            alternative = alternative,
            n = n, tol = tol
        )
    })
    number_of_obj <- length(row_ids)
    matrix_counts <- matrix(
        overlaps["count", , drop = TRUE],
        ncol = number_of_obj
    )
    matrix_pvals <- matrix(
        overlaps["pvalue", , drop = TRUE],
        ncol = number_of_obj
    )
    matrix_signs <- matrix(
        overlaps["sign", , drop = TRUE],
        ncol = number_of_obj
    )
    list(
        counts = matrix_counts,
        pvalue = matrix_pvals,
        signs = matrix_signs
    )
}

#' Rank-Rank Hypergeometric Overlap Test
#'
#' This function just extract the significant items based on RRHO analysis
#' results.
#' @param rrho_obj a object returned by [run_rrho()]
#' @param quadrant one of c("up-up", "down-down") or both, controls which kind
#' of iterms should be extracted.
#' @return a list
#' @export
rrho_sig_terms <- function(rrho_obj, quadrant = c("up-up", "down-down")) {
    res <- vector("list", length(quadrant))
    hypermat_signed <- rrho_obj$hypermat * rrho_obj$hypermat_signs
    hypermat_index <- seq(
        1L, rrho_obj$nitems,
        by = rrho_obj$stepsize
    )
    res <- lapply(quadrant, function(x) {
        quadrant_dir <- strsplit(x, "-")[[1L]]
        # transform logical index to integer value and omit NA value
        quadrant_row_index <- which(
            get_direction(rrho_obj$ranked_list$list1[hypermat_index]) ==
                quadrant_dir[[1L]]
        )
        quadrant_col_index <- which(
            get_direction(rrho_obj$ranked_list$list2[hypermat_index]) ==
                quadrant_dir[[2L]]
        )
        quadrant_hypermat <- hypermat_signed[
            quadrant_row_index, quadrant_col_index
        ]
        if (!length(quadrant_hypermat)) {
            return(NULL)
        }
        quadrant_max_value <- max(quadrant_hypermat, na.rm = TRUE)
        quadrant_sig_coord <- which(
            abs(quadrant_hypermat - quadrant_max_value) <
                sqrt(.Machine$double.eps),
            arr.ind = TRUE
        )
        sig_coord <- integer(2L)
        sig_coord[[1L]] <- quadrant_row_index[quadrant_sig_coord[1L, "row"]]
        sig_coord[[2L]] <- quadrant_col_index[quadrant_sig_coord[1L, "col"]]
        row_index <- switch(quadrant_dir[[1L]],
            up = seq_len(hypermat_index[[sig_coord[[1L]]]]),
            down = seq_len(rrho_obj$nitems) >=
                hypermat_index[[sig_coord[[1L]]]]
        )
        col_index <- switch(quadrant_dir[[2L]],
            up = seq_len(hypermat_index[[sig_coord[[2L]]]]),
            down = seq_len(rrho_obj$nitems) >=
                hypermat_index[[sig_coord[[2L]]]]
        )
        list(
            list1 = names(rrho_obj$ranked_list$list1)[row_index],
            list2 = names(rrho_obj$ranked_list$list2)[col_index],
            items = intersect(
                names(rrho_obj$ranked_list$list1)[row_index],
                names(rrho_obj$ranked_list$list2)[col_index]
            ),
            pvalue = rrho_obj$hypermat_pvalue[
                sig_coord[[1L]], sig_coord[[2L]]
            ]
        )
    })
    names(res) <- quadrant
    res
}
#' Rank-Rank Hypergeometric Overlap Map heatmap
#'
#' @param rrho_obj a object returned by [run_rrho()]
#' @param labels Not used currently.
#' @param col A vector of colors if the color mapping is discrete or a color
#' mapping function if the matrix is continuous numbers (should be generated by
#' `circlize::colorRamp2`). If the matrix is continuous, the value can also be a
#' vector of colors so that colors can be interpolated. Pass to `ColorMapping`.
#' For more details and examples, please refer to [ComplexHeatmap::Heatmap]
#' @param ... other parameters passed to [ComplexHeatmap::Heatmap]
#' @export
rrho_heatmap <- function(rrho_obj, labels, col = NULL, ...) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        rlang::abort(
            "ComplexHeatmap must be installed to use this function."
        )
    }
    if (!requireNamespace("circlize", quietly = TRUE)) {
        rlang::abort(
            "circlize must be installed to use this function."
        )
    }
    hypermat_signed <- t(rrho_obj$hypermat * rrho_obj$hypermat_signs)
    hypermat_signed <- hypermat_signed[
        rev(seq_len(ncol(hypermat_signed))), ,
        drop = FALSE
    ]
    hypermat_index <- seq(
        1L, rrho_obj$nitems,
        by = rrho_obj$stepsize
    )

    # Since now the matrix is transposed
    # row corresponding to list2
    # column corresponding to list1

    # row ranked list - left bar
    row_ranked_list <- rev(rrho_obj$ranked_list$list2[hypermat_index])
    row_label <- character(length(row_ranked_list))
    if (any(row_ranked_list > 0L)) {
        row_label[[length(row_ranked_list)]] <- "up"
    }
    if (any(row_ranked_list < 0L)) {
        row_label[[1L]] <- "down"
    }
    row_label_graphic <- list(
        up = function(x, y, w, h) {
            grid::grid.text(
                "up-regulated", x, y,
                hjust = 1L, rot = -90L,
                grid::gpar(fontface = "bold")
            )
        },
        down = function(x, y, w, h) {
            grid::grid.text(
                "down-regulated", x, y,
                hjust = 0L, rot = -90L,
                grid::gpar(fontface = "bold")
            )
        }
    )

    # column ranked list - bottom bar
    column_ranked_list <- rrho_obj$ranked_list$list1[hypermat_index]
    column_label <- character(length(column_ranked_list))
    if (any(column_ranked_list > 0L)) {
        column_label[[1L]] <- "up"
    }
    if (any(column_ranked_list < 0L)) {
        column_label[[length(column_ranked_list)]] <- "down"
    }
    column_label_graphic <- list(
        up = function(x, y, w, h) {
            grid::grid.text(
                "up-regulated", x, y,
                hjust = 0L, grid::gpar(fontface = "bold")
            )
        },
        down = function(x, y, w, h) {
            grid::grid.text(
                "down-regulated", x, y,
                hjust = 1L, grid::gpar(fontface = "bold")
            )
        }
    )

    # split parameters
    row_split <- factor(
        get_direction(row_ranked_list),
        levels = c("down", "up")
    )
    column_split <- factor(
        get_direction(column_ranked_list),
        levels = c("up", "down")
    )

    if (is.null(col)) {
        col <- circlize::colorRamp2(
            seq(
                min(hypermat_signed, na.rm = TRUE),
                max(hypermat_signed, na.rm = TRUE),
                length.out = 9L
            ),
            c(
                "#00007F", "blue", "#007FFF", "cyan",
                "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"
            )
        )
    }
    if (abs(rrho_obj$log_base - exp(1L)) < sqrt(.Machine$double.eps)) {
        legend_name <- ""
    } else {
        legend_name <- rrho_obj$log_base
    }
    ComplexHeatmap::Heatmap(
        hypermat_signed,
        col = col,
        name = paste0("-log", legend_name, "(P-value)"),
        column_title = "Rank Rank Hypergeometric Overlap Map",
        column_title_gp = grid::gpar(fontface = "bold"),
        row_title = NULL,
        left_annotation = ComplexHeatmap::rowAnnotation(
            list2_label = ComplexHeatmap::anno_customize(
                row_label,
                graphics = row_label_graphic,
                verbose = FALSE, border = FALSE
            ),
            list2 = ComplexHeatmap::anno_simple(
                row_ranked_list,
                col = circlize::colorRamp2(
                    c(-1L, 0L, 1L) * max(abs(row_ranked_list), na.rm = TRUE),
                    c("blue", "#EEEEEE", "red")
                )
            ),
            show_annotation_name = FALSE,
            show_legend = FALSE
        ),
        bottom_annotation = ComplexHeatmap::columnAnnotation(
            list1 = ComplexHeatmap::anno_simple(
                column_ranked_list,
                col = circlize::colorRamp2(
                    c(-1L, 0L, 1L) * max(
                        abs(column_ranked_list),
                        na.rm = TRUE
                    ),
                    c("blue", "#EEEEEE", "red")
                )
            ),
            list1_label = ComplexHeatmap::anno_customize(
                column_label,
                graphics = column_label_graphic,
                verbose = FALSE, border = FALSE
            ),
            show_annotation_name = FALSE,
            show_legend = FALSE
        ),
        row_split = row_split,
        cluster_rows = FALSE,
        cluster_row_slices = FALSE,
        column_split = column_split,
        cluster_columns = FALSE,
        cluster_column_slices = FALSE,
        heatmap_legend_param = list(
            title_position = "leftcenter-rot"
        ),
        ...
    )
}

#' Multiple hypothesis corrections for RRHO analysis
#'
#' Computes the significance of the agreements between lists as returned by RRHO
#' using resampling or by Benjamini-Yekutieli correction.
#'
#' @param rrho_obj a object returned by [run_rrho()]
#' @param method one of `c("BY", "permutation")` indicates which method to use
#' to correct multiple hypothesis
#' @param perm if use "permutation", how many times it uses.
#' @seealso
#' <https://academic.oup.com/nar/article/38/17/e169/1033168#82642617>
#' @export
rrho_corrected_pval <- function(rrho_obj, method = NULL, perm = 200L) {
    method <- match.arg(method, c("BY", "permutation"))
    if (identical(method, "BY")) {
        ## Convert hypermat to a vector and Benjamini Yekutieli FDR correct
        hypermatvec <- matrix(
            rrho_obj$hypermat_pvalue,
            nrow = nrow(rrho_obj$hypermat_pvalue) *
                ncol(rrho_obj$hypermat_pvalue),
            ncol = 1L
        )
        hypermat_byvec <- stats::p.adjust(hypermatvec, method = "BY")
        hyperby_mat <- matrix(
            hypermat_byvec,
            nrow = nrow(rrho_obj$hypermat_pvalue),
            ncol = ncol(rrho_obj$hypermat_pvalue)
        )
        pvalue <- min(hyperby_mat, na.rm = TRUE)
        list(
            pvalue = pvalue,
            log_pvalue = -log(pvalue, base = rrho_obj$log_base)
        )
    } else {
        perm_pvalues <- future.apply::future_vapply(
            seq_len(perm), function(i) {
                perm_rrho(rrho_obj)
            },
            FUN.VALUE = numeric(1L),
            USE.NAMES = FALSE,
            future.globals = list(rrho_obj = rrho_obj),
            future.seed = TRUE
        )
        pecdf <- function(x) {
            min(stats::ecdf(perm_pvalues)(x) + 1L / perm, 1L)
        }
        pecdf(min(rrho_obj$hypermat_pvalue, na.rm = TRUE))
    }
}

perm_items <- function(list) {
    perm_index <- sample.int(length(list), replace = FALSE)
    names(list) <- names(list)[perm_index]
    list
}
perm_rrho <- function(rrho_obj) {
    list1 <- perm_items(rrho_obj$ranked_list$list1)
    list2 <- perm_items(rrho_obj$ranked_list$list2)
    res <- calculate_hyper_overlap(
        list1, list2,
        n = length(list1),
        stepsize = rrho_obj$stepsize,
        alternative = rrho_obj$alternative
    )
    min(res$pvalue, na.rm = TRUE)
}

get_direction <- function(x) {
    data.table::fcase(
        x < 0L, "down",
        x > 0L, "up"
    )
}
