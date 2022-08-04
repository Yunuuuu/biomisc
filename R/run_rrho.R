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
#     sample1, sample2, n, 10
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
# # set.seed(1)
# rrho_correct_pval(res3, "perm", 10)
# rrho_correct_pval(res3, "perm", 10)
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
#     -log(res3$hyper_pvalue),
#     res4$hypermat
# ))
# all(dplyr::near(
#     res3$hyper_counts,
#     res4$hypermat.counts
# ))
# all(dplyr::near(
#     sign(res3$hyper_metric),
#     res4$hypermat.signs
# ))
# upgenes1 <- rrho_sig_terms(res3)
# jet.colors <- colorRampPalette(
#     c(
#         "#00007F", "blue", "#007FFF", "cyan",
#         "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"
#     )
# )
# layout(matrix(c(rep(1, 5), 2), 1, 6, byrow = TRUE))
# image(res3$hyper_metric,
#     xlab = "", ylab = "", col = jet.colors(100),
#     axes = FALSE, main = "Rank Rank Hypergeometric Overlap Map"
# )
# image(-log(res3$hyper_pvalue) * res4$hypermat.signs,
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
#' @param log_base Normally, `hyper_metric` in the results are
#' transformed by logarithm, this control the logarithm base. Just like the
#' `base` parameter in [base::log] function. Default: `10L`.
#' @references
#' <https://systems.crump.ucla.edu/rankrank/PlaisierSupplemetaryData-SupplementaryMethods_UsersGuide.pdf>
#' @export
run_rrho <- function(list1, list2, stepsize = NULL, log_base = 10L) {
    rrho_data <- set_rrho_list(list1, list2)
    if (is.null(stepsize)) {
        stepsize <- ceiling(sqrt(length(rrho_data$list2)))
    }
    ## DO Rank Rank Hypergeometric Overlap
    hyper_res <- calculate_hyper_overlap(
        rrho_data$list1, rrho_data$list2,
        n = length(rrho_data$list1),
        stepsize = stepsize
    )
    hyper_metric <- -log(hyper_res$pvalue, base = log_base) *
        hyper_res$signs

    list(
        hyper_metric = hyper_metric,
        hyper_pvalue = hyper_res$pvalue,
        hyper_counts = hyper_res$counts,
        rrho_data = rrho_data,
        stepsize = stepsize,
        log_base = log_base
    )
}
set_rrho_list <- function(list1, list2) {
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
    list(
        list1 = sort(list1, decreasing = TRUE),
        list2 = sort(list2, decreasing = TRUE)
    )
}

# Manuscript definition
# N: the length of list1 - population size
# M: the current rank step pair of sample1 - target population
# s: the current rank step pair of sample2 - sample size
# k: overlapping length

# `phyper` function definition
# count = q = k: overlapping length
# m = M: the current rank step pair of sample1 - target population
# n = N-M: non-target population
# k = s: the current rank step pair of sample2 - sample size

## Compute the overlaps between two *numeric* lists:
hyper_test <- function(sample1, sample2, n) {
    count <- length(intersect(sample1, sample2))
    m <- length(sample1)
    k <- length(sample2)
    pvalue <- stats::phyper(
        q = count - 1L, m = m, n = n - m + 1L,
        k = k, lower.tail = FALSE, log.p = FALSE
    )
    sign <- sign(count - m * k / n)
    c(
        count = count,
        pvalue = pvalue,
        sign = sign
    )
}

calculate_hyper_overlap <- function(list1, list2, n, stepsize) {
    row_ids <- seq(1L, length(list1), by = stepsize)
    col_ids <- seq(1L, length(list2), by = stepsize)
    indexes <- expand.grid(
        row_ids = row_ids,
        col_ids = col_ids
    )
    indexes <- as.matrix(indexes)
    overlaps <- apply(indexes, 1L, function(x) {
        hyper_test(
            names(list1)[seq_len(x[["row_ids"]])],
            names(list2)[seq_len(x[["col_ids"]])],
            n = n
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
#' @param quadrant one or more items in c("up-up", "down-down", "up-down",
#' "down-up"), controls which kind of iterms should be extracted.
#' @return a list
#' @export
rrho_sig_terms <- function(rrho_obj, quadrant = c("up-up", "down-down")) {
    stopifnot(
        all(quadrant %in% c("up-up", "down-down", "up-down", "down-up"))
    )
    res <- vector("list", length(quadrant))
    rrho_list1_index <- seq(
        1L, length(rrho_obj$rrho_data$list1),
        by = rrho_obj$stepsize
    )
    rrho_list2_index <- seq(
        1L, length(rrho_obj$rrho_data$list2),
        by = rrho_obj$stepsize
    )
    res <- lapply(quadrant, function(x) {
        quadrant_dir <- strsplit(x, "-")[[1L]]

        # integer index relative to rrho_obj$hyper_metric and rrho_obj$hyper_pvalue
        # the row index is relative to `rrho_list1_index`
        # the column index is relative to `rrho_list2_index`
        quadrant_row_index <- which(
            get_direction(rrho_obj$rrho_data$list1[rrho_list1_index]) ==
                quadrant_dir[[1L]]
        )
        quadrant_col_index <- which(
            get_direction(rrho_obj$rrho_data$list2[rrho_list2_index]) ==
                quadrant_dir[[2L]]
        )
        quadrant_hyper_metric <- rrho_obj$hyper_metric[
            quadrant_row_index, quadrant_col_index
        ]
        quadrant_hyper_counts <- rrho_obj$hyper_counts[
            quadrant_row_index, quadrant_col_index
        ]
        if (!length(quadrant_hyper_metric)) {
            return(NULL)
        }

        # For `quadrant_sig_coord`
        # integer index relative to `quadrant_hyper_metric`
        # the row index is relative to `quadrant_row_index`
        # the column index is relative to `quadrant_col_index`
        if (x %in% c("up-up", "down-down")) {
            quadrant_sig_value <- max(quadrant_hyper_metric, na.rm = TRUE)
            if (quadrant_sig_value <= 0L) {
                return(NULL)
            }
            if (is.infinite(quadrant_sig_value)) {
                quadrant_sig_coord <- which(
                    is.infinite(quadrant_hyper_metric) &
                        sign(quadrant_hyper_metric) > 0L,
                    arr.ind = TRUE
                )
            } else {
                quadrant_sig_coord <- which(
                    abs(quadrant_hyper_metric - quadrant_sig_value) <
                        sqrt(.Machine$double.eps),
                    arr.ind = TRUE
                )
            }
        } else if (x %in% c("up-down", "down-up")) {
            quadrant_sig_value <- min(quadrant_hyper_metric, na.rm = TRUE)
            if (quadrant_sig_value >= 0L) {
                return(NULL)
            }
            if (is.infinite(quadrant_sig_value)) {
                quadrant_sig_coord <- which(
                    is.infinite(quadrant_hyper_metric) &
                        sign(quadrant_hyper_metric) < 0L,
                    arr.ind = TRUE
                )
            } else {
                quadrant_sig_coord <- which(
                    abs(quadrant_hyper_metric - quadrant_sig_value) <
                        sqrt(.Machine$double.eps),
                    arr.ind = TRUE
                )
            }
        }
        quadrant_sig_coord <- quadrant_sig_coord[
            which.max(quadrant_hyper_counts[quadrant_sig_coord]), ,
            drop = TRUE
        ]
        quadrant_sig_coord <- c(
            switch(quadrant_dir[[1L]],
                up = max(quadrant_sig_coord[[1L]]),
                down = min(quadrant_sig_coord[[1L]])
            ),
            switch(quadrant_dir[[2L]],
                up = max(quadrant_sig_coord[[2L]]),
                down = min(quadrant_sig_coord[[2L]])
            )
        )
        sig_coord <- c(
            rrho_list1_index[quadrant_row_index[
                quadrant_sig_coord[[1L]]
            ]],
            rrho_list2_index[quadrant_col_index[
                quadrant_sig_coord[[2L]]
            ]]
        )
        list1_index <- switch(quadrant_dir[[1L]],
            up = seq_len(sig_coord[[1L]]),
            down = seq_along(rrho_obj$rrho_data$list1) >= sig_coord[[1L]]
        )
        list2_index <- switch(quadrant_dir[[2L]],
            up = seq_len(sig_coord[[2L]]),
            down = seq_along(rrho_obj$rrho_data$list2) >= sig_coord[[2L]]
        )
        list(
            sig_item1 = names(rrho_obj$rrho_data$list1)[list1_index],
            sig_item2 = names(rrho_obj$rrho_data$list2)[list2_index],
            common_items = intersect(
                names(rrho_obj$rrho_data$list1)[list1_index],
                names(rrho_obj$rrho_data$list2)[list2_index]
            ),
            pvalue = rrho_obj$hyper_pvalue[
                quadrant_row_index[quadrant_sig_coord[[1L]]],
                quadrant_col_index[quadrant_sig_coord[[2L]]]
            ],
            hyper_metric = rrho_obj$hyper_metric[
                quadrant_row_index[quadrant_sig_coord[[1L]]],
                quadrant_col_index[quadrant_sig_coord[[2L]]]
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
    # heat_matrix <- t(rrho_obj$hyper_metric)
    heat_matrix <- t(-log(rrho_obj$hyper_pvalue, base = rrho_obj$log_base))
    heat_matrix <- heat_matrix[
        rev(seq_len(nrow(heat_matrix))), ,
        drop = FALSE
    ]

    rrho_list1_index <- seq(
        1L, length(rrho_obj$rrho_data$list1),
        by = rrho_obj$stepsize
    )
    rrho_list2_index <- seq(
        1L, length(rrho_obj$rrho_data$list2),
        by = rrho_obj$stepsize
    )

    # Since now the matrix is transposed
    # row corresponding to list2
    # column corresponding to list1

    # row ranked list - left bar
    row_ranked_list <- rev(rrho_obj$rrho_data$list2[rrho_list2_index])
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
    column_ranked_list <- rrho_obj$rrho_data$list1[rrho_list1_index]
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
                min(heat_matrix[is.finite(heat_matrix)], na.rm = TRUE),
                max(heat_matrix[is.finite(heat_matrix)], na.rm = TRUE),
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
    legend_name <- paste0("-log", legend_name, "(P-value)")
    ComplexHeatmap::Heatmap(
        heat_matrix,
        col = col, name = legend_name,
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
#' @param perm if use "permutation", how many permutation times it uses.
#' @seealso
#' <https://academic.oup.com/nar/article/38/17/e169/1033168#82642617>
#' @export
rrho_correct_pval <- function(rrho_obj, method = NULL, perm = 200L) {
    method <- match.arg(method, c("BY", "permutation"))
    if (identical(method, "BY")) {
        ## Convert hypermat to a vector and Benjamini Yekutieli FDR correct
        hyper_pvalue_by <- stats::p.adjust(c(rrho_obj$hyper_pvalue), method = "BY")
        hyper_pvalue_by <- matrix(
            hyper_pvalue_by,
            nrow = nrow(rrho_obj$hyper_pvalue),
            ncol = ncol(rrho_obj$hyper_pvalue)
        )
        list(
            hyper_metric_by = -log(hyper_pvalue_by, base = rrho_obj$log_base),
            hyper_pvalue_by = hyper_pvalue_by
        )
    } else {
        perm_hyper_metric <- future.apply::future_vapply(
            seq_len(perm), function(i) {
                perm_rrho(rrho_obj)
            },
            FUN.VALUE = numeric(1L),
            USE.NAMES = FALSE,
            future.globals = list(rrho_obj = rrho_obj),
            future.seed = TRUE
        )
        pecdf <- stats::ecdf(perm_hyper_metric)
        list(
            ecdf = pecdf,
            pvalue_perm = 1 - min(
                pecdf(max(rrho_obj$hyper_metric, na.rm = TRUE)) + 1L / perm,
                1L
            )
        )
    }
}
perm_items <- function(list) {
    perm_index <- sample.int(length(list), replace = FALSE)
    names(list) <- names(list)[perm_index]
    list
}
perm_rrho <- function(rrho_obj) {
    list1 <- perm_items(rrho_obj$rrho_data$list1)
    list2 <- perm_items(rrho_obj$rrho_data$list2)
    hyper_res <- calculate_hyper_overlap(
        list1, list2,
        n = length(list1),
        stepsize = rrho_obj$stepsize
    )
    hyper_metric <- -log(hyper_res$pvalue, base = rrho_obj$log_base) *
        hyper_res$signs
    max(hyper_metric, na.rm = TRUE)
}

get_direction <- function(x) {
    data.table::fcase(
        x < 0L, "down",
        x > 0L, "up"
    )
}
