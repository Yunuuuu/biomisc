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
# res1 <- rrho_hyper_overlap(
#     names(sample1), names(sample2), length(sample1), 1L
# )
# res2 <- RRHO:::numericListOverlap(
#     names(sample1), names(sample2), 1,
#     alternative = "enrichment"
# )
# all(dplyr::near(
#     res1$counts, res2$counts
# ))
# all(dplyr::near(
#     -log(res1$pvalue), res2$log.pval
# ))
#
# rrho_res <- run_rrho(sample1, sample2, 1)
# rrho_heatmap(rrho_res)
# set.seed(1)
# rrho_correct_pval(rrho_res, "perm", 10)
# rrho_correct_pval(rrho_res, "perm", 10)
# set.seed(1)
# rrho_correct_pval(rrho_res, "perm", 10)
# rrho_correct_pval(rrho_res, "perm", 10)
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
# upgenes1 <- rrho_sig_items(res3)
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
#' The function tests for significant overlap between two sorted lists using the
#' method in the reference. Modified from <https://github.com/johnros/RRHO> by
#' making adjustment with official RRHO described in
#' <https://academic.oup.com/nar/article/38/17/e169/1033168#82642652>.
#'
#' @param list1,list2 A named numeric vector, For differential gene expression,
#' values are often `-log10(P-value) * sign(effect)`. `list1` will be regarded
#' as the reference populations.
#' @param stepsize Controls the resolution of the test: how many items between
#' any two overlap tests.
#' @param correction A string indicates the correction methods when `list1` and
#' `list2` have different length. One of `c("common", "length")`, if "common",
#' the `scale_size` will be equal to 1L.
#' @param log_base Normally, `hyper_metric` in the results are transformed by
#' logarithm, this control the logarithm base. Just like the `base` parameter in
#' [base::log] function. Default: `10L`.
#' @return A list of class "rrho" with the following fields:
#' \describe{
#'  \item{hyper_metric}{the metric of Rank-Rank Hypergeometric Overlap analysis,
#'  which is the basis of the RRHO heatmap and permutation test, usually equals
#'  to `abs(log(hyper_pvalue)) * hyper_signs`.}
#'  \item{hyper_pvalue}{The Hypergeometric test cumulative pvalue.}
#'  \item{hyper_signs}{the signs of `hyper_metric`, negative values indicate
#'  under-enrichment, and positive values indicate over-enrichment.}
#'  \item{hyper_counts}{the overlapping counts number in each test.}
#'  \item{rrho_data}{the data used in RRHO analysis: both `list1` (all items or
#'  only keep items also in `list2` base on "correction" argument) and `list2`
#'  (will only keep items in `list1`) and sort them based on their value.}
#'  \item{stepsize}{the provided `stepsize` parameter.}
#'  \item{log_base}{the provided `log_base` parameter.}
#' }
#' @examples
#' n <- 200
#' sample1 <- rnorm(n)
#' sample2 <- rnorm(n)
#' names(sample1) <- names(sample2) <- paste0("gene", seq_len(n))
#' rrho_res <- biomisc::run_rrho(sample1, sample2, 1)
#' @references
#' * <https://academic.oup.com/nar/article/38/17/e169/1033168#82642652>
#' * <https://systems.crump.ucla.edu/rankrank/PlaisierSupplemetaryData-SupplementaryMethods_UsersGuide.pdf>
#' @export
#' @rdname run_rrho
run_rrho <- function(list1, list2, stepsize = NULL, correction = NULL, log_base = 10L) {
    correction <- match.arg(correction, c("common", "length"))
    assert_class(stepsize, is_scalar_numeric, "scalar numeric", null_ok = TRUE)
    rrho_data <- set_rrho_list(list1, list2, correction = correction)
    if (is.null(stepsize)) {
        stepsize <- as.integer(sqrt(min(lengths(rrho_data)[1:2])))
    } else {
        stepsize <- max(1L, min(as.integer(stepsize), lengths(rrho_data)[1:2]))
    }
    ## DO Rank Rank Hypergeometric Overlap
    hyper_res <- rrho_hyper_overlap(
        names(rrho_data$list1),
        names(rrho_data$list2),
        stepsize = stepsize
    )
    new_rrho(list(
        hyper_pvalue = hyper_res$pvalue,
        hyper_signs = hyper_res$signs,
        hyper_counts = hyper_res$counts,
        rrho_data = rrho_data,
        stepsize = stepsize,
        log_base = log_base
    ))
}

new_rrho <- function(list) {
    list$hyper_metric <- rrho_metric(
        list$hyper_pvalue,
        list$hyper_signs,
        list$rrho_data$scale_size,
        list$log_base
    )
    structure(list, class = "rrho")
}

rrho_metric <- function(pvalue, signs, scale_size, log_base) {
    signs * scale_size * abs(log(pvalue, log_base))
}

#' @param x An object returned by [run_rrho()]
#' @param ... Not used currently
#' @export
#' @rdname run_rrho
print.rrho <- function(x, ...) {
    cat(strwrap(
        "Rank-Rank Hypergeometric Overlap analysis",
        indent = 0, exdent = 2
    ), sep = "\n")
    sig_spots <- rrho_sig_spot_internal(x)
    cat(strwrap(
        paste(
            "The maximal RRHO metrix",
            paste(sprintf(
                "%.2g (%s-%s)", sig_spots$hyper_metric,
                sig_spots$list1, sig_spots$list2
            ), collapse = ", "),
            sep = ": "
        ),
        indent = 2, exdent = 2
    ), sep = "\n")
    if (is.integer(x$rrho_data$scale_size)) {
        sprintf_term <- "%d"
    } else {
        sprintf_term <- "%.2f"
    }
    cat(strwrap(
        sprintf(
            paste0("RRHO metrix scale_size: ", sprintf_term),
            x$rrho_data$scale_size
        ),
        indent = 2, exdent = 2
    ), sep = "\n")
    cat(strwrap(
        sprintf("Analysis with stepsize: %d", x$stepsize),
        indent = 2, exdent = 2
    ), sep = "\n")
}

rrho_sig_spot_internal <- function(rrho_obj) {
    rrho_list1_index <- seq.int(
        rrho_obj$stepsize, length(rrho_obj$rrho_data$list1),
        by = rrho_obj$stepsize
    )
    rrho_list2_index <- seq.int(
        rrho_obj$stepsize, length(rrho_obj$rrho_data$list2),
        by = rrho_obj$stepsize
    )
    idx <- which(
        rrho_obj$hyper_pvalue == min(rrho_obj$hyper_pvalue, na.rm = TRUE),
        arr.ind = TRUE
    )
    out <- data.table::data.table(
        hyper_metric = rrho_obj$hyper_metric[idx],
        hyper_pvalue = rrho_obj$hyper_pvalue[idx],
        list1 = rrho_get_direction(rrho_obj$rrho_data$list1[rrho_list1_index])[
            idx[, 1L, drop = TRUE]
        ],
        list2 = rrho_get_direction(rrho_obj$rrho_data$list2[rrho_list2_index])[
            idx[, 2L, drop = TRUE]
        ]
    )
    unique(out)
}

set_rrho_list <- function(list1, list2, correction) {
    # remove NA value
    list1 <- list1[!is.na(list1)]
    list2 <- list2[!is.na(list2)]

    # keep items in both lists
    common_names <- intersect(names(list1), names(list2))
    list2_filtered <- list2[common_names]
    if (identical(correction, "common")) {
        list1 <- list1[common_names]
        scale_size <- 1L
        cli::cli_inform(
            "Finding {length(common_names)} genes shared by {.field list1} and {.field list2}"
        )
    } else if (identical(correction, "length")) {
        scale_size <- length(list2_filtered) / length(list1)
        cli::cli_inform(
            "Removing {length(list2) - length(list2_filtered)} genes from {.field list2} not in {.field list1}"
        )
    }

    list(
        list1 = sort(list1, decreasing = TRUE),
        list2 = sort(list2_filtered, decreasing = TRUE),
        scale_size = scale_size
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

## Compute the overlaps between two *character* atomic vector:
hyper_test <- function(sample1, sample2, n) {
    count <- length(intersect(sample1, sample2))
    m <- length(sample1)
    k <- length(sample2)

    # under-enrichment
    if (count <= m / n * k) { # fix error: NAs produced by integer overflow
        sign <- -1L
        pvalue <- stats::phyper(
            q = count, m = m, n = n - m,
            k = k, lower.tail = TRUE, log.p = FALSE
        )
    } else {
        # over-enrichment
        # since lower.tail = FALSE won't include point estimation in `count`
        # value, so we just subtract one to include the point estimation
        # Also since count > m * k / n, `count` won't be able to equal to zero
        # it's safe to just subtract one
        sign <- 1L
        pvalue <- stats::phyper(
            q = count - 1L, m = m, n = n - m,
            k = k, lower.tail = FALSE, log.p = FALSE
        )
    }

    c(count, pvalue, sign)
}

rrho_hyper_overlap <- function(sample1, sample2, stepsize) {
    n <- length(sample1)
    row_ids <- seq.int(stepsize, length(sample1), by = stepsize)
    col_ids <- seq.int(stepsize, length(sample2), by = stepsize)
    indexes <- expand.grid(
        row_ids = row_ids,
        col_ids = col_ids
    )
    row_ind <- indexes[["row_ids"]]
    col_ind <- indexes[["col_ids"]]
    p <- progressr::progressor(steps = nrow(indexes) / 100L)
    overlaps <- future.apply::future_lapply(
        seq_len(nrow(indexes)),
        function(i) {
            if (i %% 100L == 0L) {
                p(message = "hyper-geometric testing")
            }
            hyper_test(
                sample1[seq_len(row_ind[[i]])],
                sample2[seq_len(col_ind[[i]])],
                n = n
            )
        }
    )
    overlaps <- data.table::transpose(overlaps)
    number_of_obj <- length(row_ids)
    matrix_counts <- matrix(
        as.integer(overlaps[[1L]]),
        nrow = number_of_obj
    )
    # check Pvalue range from [0L, 1L]
    # if (any(overlaps[[2L]] > 1L | overlaps[[2L]] < 0L)) {
    #     cli::cli_abort("Something wrong when calculating Hypergeometric Distribution pvalue")
    # }
    matrix_pvals <- matrix(
        overlaps[[2L]],
        nrow = number_of_obj
    )
    matrix_signs <- matrix(
        as.integer(overlaps[[3L]]),
        nrow = number_of_obj
    )
    list(
        counts = matrix_counts,
        pvalue = matrix_pvals,
        signs = matrix_signs
    )
}

#' Rank-Rank Hypergeometric Overlap significant items
#'
#' This function just extract the significant items based on RRHO analysis
#' results.
#' @param rrho_obj An object returned by [run_rrho()]
#' @param quadrant one or more items in `c("up-up", "down-down", "up-down",
#' "down-up")`, controls which quadrant of items should be extracted.
#' @details The highest intensity point on the resulting Rank-Rank
#' Hypergeometric Overlap map corresponds to the pair of rank thresholds where
#' the observed statistical overlap between the two gene-expression profiles is
#' the strongest statistically. In other words, the highest intensity point
#' represents the optimal overlap between the profiles in that this is the
#' overlap that is least likely to occur by chance. The coordinates of the
#' highest intensity point are the rank in each experiment above which are the
#' most statistically significant set of overlapping genes.
#' @return a list of object for each quadrant with class "rrho_sig" which
#' contains the following fields:
#' \describe{
#'  \item{sig_item1}{the sinificant items in list1.}
#'  \item{sig_item2}{the sinificant items in list2.}
#'  \item{common_items}{the intersected items between `sig_item1` and
#'  `sig_item2`.}
#'  \item{hyper_metric}{the corresponding RRHO metric in the significant
#'  coordinate.}
#'  \item{hyper_pvalue}{the corresponding RRHO pvalue in in the significant
#'  coordinate.}
#' }
#' @references
#' <https://academic.oup.com/nar/article/38/17/e169/1033168#82642652>
#' @seealso [rrho_sig_spot]
#' @examples
#' n <- 200
#' sample1 <- rnorm(n)
#' sample2 <- rnorm(n)
#' names(sample1) <- names(sample2) <- paste0("gene", seq_len(n))
#' rrho_res <- biomisc::run_rrho(sample1, sample2, 1)
#' biomisc::rrho_sig_items(rrho_res)
#' @export
#' @rdname rrho_sig_items
rrho_sig_items <- function(rrho_obj, quadrant = c("up-up", "down-down")) {
    assert_rrho(rrho_obj)
    stopifnot(
        all(quadrant %chin% c("up-up", "down-down", "up-down", "down-up"))
    )
    rrho_list1_index <- seq.int(
        rrho_obj$stepsize, length(rrho_obj$rrho_data$list1),
        by = rrho_obj$stepsize
    )
    rrho_list2_index <- seq.int(
        rrho_obj$stepsize, length(rrho_obj$rrho_data$list2),
        by = rrho_obj$stepsize
    )
    res <- lapply(quadrant, function(x) {
        quadrant_dir <- strsplit(x, "-")[[1L]]

        # integer index relative to rrho_obj$hyper_metric and
        # rrho_obj$hyper_pvalue, namely the same length with corresponding dim.
        # the row index is relative to `rrho_list1_index`
        # the column index is relative to `rrho_list2_index`
        quadrant_row_index <- which(
            rrho_get_direction(rrho_obj$rrho_data$list1[rrho_list1_index]) ==
                quadrant_dir[[1L]]
        )
        quadrant_col_index <- which(
            rrho_get_direction(rrho_obj$rrho_data$list2[rrho_list2_index]) ==
                quadrant_dir[[2L]]
        )
        quadrant_hyper_metric <- rrho_obj$hyper_metric[
            quadrant_row_index, quadrant_col_index
        ]
        quadrant_hyper_counts <- rrho_obj$hyper_counts[
            quadrant_row_index, quadrant_col_index
        ]
        quadrant_hyper_pvalue <- rrho_obj$hyper_pvalue[
            quadrant_row_index, quadrant_col_index
        ]
        if (!length(quadrant_hyper_pvalue)) {
            return(NULL)
        }

        # For `quadrant_sig_coord`
        # integer index relative to `quadrant_hyper_pvalue`
        # the row index is relative to `quadrant_row_index`
        # the column index is relative to `quadrant_col_index`

        # for "up-down" and "down-up" quadrant, under-enrichment means
        # over-enrichment, so we should find the minimal value and ensue it
        # is negative
        quadrant_sig_coord <- which(
            quadrant_hyper_pvalue == min(quadrant_hyper_pvalue, na.rm = TRUE),
            arr.ind = TRUE
        )
        quadrant_sig_sign <- quadrant_hyper_metric[quadrant_sig_coord]
        if (quadrant_dir[[1L]] == quadrant_dir[[2L]]) {
            if (quadrant_sig_sign <= 0L) {
                return(NULL)
            }
        } else {
            if (quadrant_sig_sign >= 0L) {
                return(NULL)
            }
        }

        # if there exists more than one significant values
        # Just keep the one with the largest number of items.
        # `quadrant_sig_coord` is a matrix with every row correspond to
        # every significant values
        # Here we find the row with the maximal counts: returns a integer atomic
        # vector, the first value is the sinificant value index of
        # quadrant_row_index, the second quadrant_col_index
        #
        # notes: if there are more than one significant values with the same
        # number of counts number, this will just keep the first one. But it's
        # not harmful for us since we usually only need the overlapping items
        quadrant_sig_coord <- quadrant_sig_coord[
            which.max(quadrant_hyper_counts[quadrant_sig_coord]), ,
            drop = TRUE
        ]
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
        sig_item1 <- names(rrho_obj$rrho_data$list1)[list1_index]
        sig_item2 <- names(rrho_obj$rrho_data$list2)[list2_index]

        structure(
            list(
                sig_item1 = sig_item1, sig_item2 = sig_item2,
                common_items = intersect(sig_item1, sig_item2),
                hyper_metric = quadrant_hyper_metric[quadrant_sig_coord],
                hyper_pvalue = quadrant_hyper_pvalue[quadrant_sig_coord]
            ),
            class = "rrho_sig"
        )
    })
    names(res) <- quadrant
    res
}

#' Rank-Rank Hypergeometric Overlap significant spot
#'
#' Extract the the significant spot and its quadrant
#'
#' @inheritParams rrho_sig_items
#' @inherit rrho_sig_items details
#' @return A data.frame, of which each row indicates the highest intensity
#'   point.
#' @seealso [rrho_sig_items]
#' @export
rrho_sig_spot <- function(rrho_obj) {
    assert_rrho(rrho_obj)
    sig_spot <- rrho_sig_spot_internal(rrho_obj)
    data.table::setDF(sig_spot)
    sig_spot
}

#' @param x An object returned by [rrho_sig_items()]
#' @param ... Not used currently
#' @export
#' @rdname rrho_sig_items
print.rrho_sig <- function(x, ...) {
    cat(strwrap(
        sprintf(
            "Finding %d significant overlapping items by list1 and list2 in RRHO analysis",
            length(x$common_items)
        ),
        indent = 0, exdent = 2
    ), sep = "\n")
    cat(strwrap(
        sprintf(
            "Significant items in list1: %d",
            length(x$sig_item1)
        ),
        indent = 2, exdent = 2
    ), sep = "\n")
    cat(strwrap(
        sprintf(
            "Significant items in list2: %d",
            length(x$sig_item2)
        ),
        indent = 2, exdent = 2
    ), sep = "\n")
    cat(strwrap(
        sprintf("The significant RRHO metrix: %.2g", x$hyper_metric),
        indent = 0, exdent = 2
    ), sep = "\n")
    cat(strwrap(
        sprintf("The significant RRHO pvalue: %.2g", x$hyper_pvalue),
        indent = 0, exdent = 2
    ), sep = "\n")
}

#' Rank-Rank Hypergeometric Overlap Map heatmap
#'
#' @inheritParams rrho_sig_items
#' @param col A vector of colors if the color mapping is discrete or a color
#' mapping function if the matrix is continuous numbers (should be generated by
#' `circlize::colorRamp2`). If the matrix is continuous, the value can also be a
#' vector of colors so that colors can be interpolated. Pass to `ColorMapping`.
#' For more details and examples, please refer to [ComplexHeatmap::Heatmap]
#' @inheritParams ComplexHeatmap::Heatmap
#' @inheritDotParams ComplexHeatmap::Heatmap -matrix -col -name -column_title
#' -column_title_gp -row_title -left_annotation -bottom_annotation -row_split
#' -cluster_rows -cluster_row_slices -column_split -cluster_columns
#' -cluster-column_slices -heatmap_legend_param -use_raster
#' @examples
#' n <- 200
#' sample1 <- rnorm(n)
#' sample2 <- rnorm(n)
#' names(sample1) <- names(sample2) <- paste0("gene", seq_len(n))
#' rrho_res <- biomisc::run_rrho(sample1, sample2, 1)
#' biomisc::rrho_heatmap(rrho_res)
#' @export
rrho_heatmap <- function(rrho_obj, col = NULL, ..., use_raster = NULL) {
    assert_rrho(rrho_obj)
    assert_pkg("ComplexHeatmap")
    assert_pkg("circlize")
    assert_class(use_raster, rlang::is_scalar_logical,
        msg = "scalar {.cls logical} value", null_ok = TRUE
    )
    heat_matrix <- t(rrho_obj$hyper_metric)
    # heat_matrix <- t(-log(rrho_obj$hyper_pvalue, base = rrho_obj$log_base))
    heat_matrix <- heat_matrix[
        rev(seq_len(nrow(heat_matrix))), ,
        drop = FALSE
    ]

    rrho_list1_index <- seq.int(
        rrho_obj$stepsize, length(rrho_obj$rrho_data$list1),
        by = rrho_obj$stepsize
    )
    rrho_list2_index <- seq.int(
        rrho_obj$stepsize, length(rrho_obj$rrho_data$list2),
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
        rrho_get_direction(row_ranked_list),
        levels = c("down", "up")
    )
    column_split <- factor(
        rrho_get_direction(column_ranked_list),
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
    if (rrho_obj$log_base == exp(1L)) {
        legend_name <- ""
    } else {
        legend_name <- rrho_obj$log_base
    }
    legend_name <- paste0("Signed |log", legend_name, "(P-value)|")
    if (is.null(use_raster)) {
        if (length(rrho_list1_index) > 500L) {
            cli::cli_inform("Setting {.code use_raster = TRUE}")
            use_raster <- TRUE
        } else {
            cli::cli_inform("Setting {.code use_raster = FALSE}")
            use_raster <- FALSE
        }
    }
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
        ...,
        use_raster = use_raster
    )
}

#' Multiple hypothesis corrections for RRHO analysis
#'
#' Computes the significance of the agreements between lists as returned by RRHO
#' using resampling or by Benjamini-Yekutieli correction.
#'
#' @inheritParams rrho_sig_items
#' @param method A string, "permutation" or other correction method (details
#' see [p.adjust()]) indicates which method to use
#' to correct multiple hypothesis
#' @param perm if use "permutation", indicates how many times should be used.
#' @param quadrant the "quadrant" to test significance, usually we want to test
#' whether the overlapping goes in the same direction (over-enrichment), which
#' means "up-up" quadrant and "down-down" quadrant. You can specify "all" to
#' test the overall significance of the RRHO map. Due to the sign convention for
#' over- or under-enrichment in RRHO design, we can test the over-enrichment of
#' "down-up" quadrant or/and "up-down" quatrant by testing the under-enrichment
#' of hotspot significance in "down-up" quadrant or/and "up-down" quatrant, the
#' Pvalue for which test whether the overlapping goes in the different
#' direction.
#' @details
#' If the Benjamini-Yekutieli corrected hypergeometric map has most significant
#' areas of absolute log10(P-value) intensity 15 or greater, then permutations
#' are likely to be significant below a permutation frequency P-value of 0.01.
#' When the corrected map maxima are less than 15, perform permutations. If the
#' two maps for comparison are made from gene lists that are of considerably
#' different lengths, one can either (a) scale the hypergeometric maps using the
#' list length correction method to account for length differences before
#' applying a Benjamini-Yekutieli correction or (b) remake the hypergeometric
#' maps using only items common to all signatures in the set of comparisons
#' making sure this restriction is not too limiting on the total gene number.
#'
#' Areas "up-up" quadrant and "down-down" quadrant correspond to overlapping
#' items that go in the same direction in both experiments. We compare the
#' maximum in these areas as a summary statistic to screen through permutation
#' RRHO maps compared to the true RRHO map. The frequency at which permutation
#' maps have a higher summary statistic than the true map is defined as the
#' permutation P-value.
#' @examples
#' n <- 200
#' sample1 <- rnorm(n)
#' sample2 <- rnorm(n)
#' names(sample1) <- names(sample2) <- paste0("gene", seq_len(n))
#' rrho_res <- biomisc::run_rrho(sample1, sample2, 1)
#' progressr::with_progress(
#'     biomisc::rrho_correct_pval(rrho_res, "permutation", 10L)
#' )
#' @seealso
#' <https://academic.oup.com/nar/article/38/17/e169/1033168#82642617>
#' <https://systems.crump.ucla.edu/rankrank/PlaisierSupplemetaryData-SupplementaryMethods_UsersGuide.pdf>
#' @export
rrho_correct_pval <- function(rrho_obj, method = "BY", perm = 200L, quadrant = c("up-up", "down-down")) {
    assert_rrho(rrho_obj)
    method <- match.arg(method, c(p.adjust.methods, "permutation"))
    if (!identical(method, "permutation")) {
        # Convert hypermat to a vector and apply Benjamini Yekutieli Pvalue
        # correction
        hyper_pvalue_by <- stats::p.adjust(
            c(rrho_obj$hyper_pvalue),
            method = method
        )
        rrho_obj$hyper_pvalue <- matrix(
            hyper_pvalue_by,
            nrow = nrow(rrho_obj$hyper_pvalue),
            ncol = ncol(rrho_obj$hyper_pvalue)
        )
        new_rrho(rrho_obj)
    } else {
        quadrant <- unique(quadrant)
        # we want to use over-enrichment for all quadrant, but we actually want
        # to a under-enrichment test in c("up-down", "down-up"), by providing a
        # `quadrant_sign`, we can easily transformed the under-enrichment into a
        # under-enrichment.
        if (all(quadrant %chin% c("up-up", "down-down"))) {
            quadrant_sign <- 1L
        } else if (all(quadrant %chin% c("up-down", "down-up"))) {
            quadrant_sign <- -1L
        } else if (!identical(quadrant, "all")) {
            cli::cli_abort(
                c(
                    "x" = "Wrong {.arg quadrant}",
                    "i" = "{.arg quadrant} should be one of {.code c(\"all\", \"up-up\", \"down-down\", \"up-down\", \"down-up\")}.", # nolint
                    "i" = "you can also specify {.code c(\"up-up\", \"down-down\")} or {.code c(\"up-down\", \"down-up\")} as a whole." # nolint
                )
            )
        }

        # define quadrant for testing
        if (!identical(quadrant, "all")) {
            rrho_list1_index <- seq.int(
                rrho_obj$stepsize, length(rrho_obj$rrho_data$list1),
                by = rrho_obj$stepsize
            )
            rrho_list2_index <- seq.int(
                rrho_obj$stepsize, length(rrho_obj$rrho_data$list2),
                by = rrho_obj$stepsize
            )
            rrho_list1_quadrant <- rrho_get_direction(
                rrho_obj$rrho_data$list1[rrho_list1_index]
            )
            rrho_list2_quadrant <- rrho_get_direction(
                rrho_obj$rrho_data$list2[rrho_list2_index]
            )
            rrho_quadrant <- outer(
                rrho_list1_quadrant, rrho_list2_quadrant,
                paste,
                sep = "-"
            )
            quadrant_idx_list <- lapply(quadrant, function(x) {
                x == rrho_quadrant
            })
        }

        # implement permutation
        p <- progressr::progressor(steps = perm)
        perm_hyper_metric <- future.apply::future_lapply(
            seq_len(perm), function(i) {
                # https://github.com/HenrikBengtsson/progressr/issues/136
                # call inner progress first then the outer progress bar
                # in case of Inner progress bars displayed after outer progress
                # finishes
                out <- rrho_hyper_overlap(
                    names(rrho_obj$rrho_data$list1)[
                        sample.int(length(rrho_obj$rrho_data$list1), replace = FALSE)
                    ],
                    names(rrho_obj$rrho_data$list2)[
                        sample.int(length(rrho_obj$rrho_data$list2), replace = FALSE)
                    ],
                    stepsize = rrho_obj$stepsize
                )
                p(message = sprintf("Permuatating %d times", i))
                rrho_metric(
                    out$pvalue, out$signs,
                    rrho_obj$rrho_data$scale_size,
                    rrho_obj$log_base
                )
            },
            future.globals = TRUE,
            future.seed = TRUE
        )
        # derive permutation summary statistics for given quadrant
        summary_stats <- vapply(perm_hyper_metric, function(hyper_metric_mat) {
            rrho_summary_stats(
                quadrant = quadrant,
                quadrant_idx_list = quadrant_idx_list,
                quadrant_sign = quadrant_sign,
                hyper_metric_mat = hyper_metric_mat
            )
        }, FUN.VALUE = numeric(1L), USE.NAMES = FALSE)
        pecdf <- stats::ecdf(summary_stats)

        actual_stats <- rrho_summary_stats(
            quadrant = quadrant,
            quadrant_idx_list = quadrant_idx_list,
            quadrant_sign = quadrant_sign,
            hyper_metric_mat = rrho_obj$hyper_metric
        )
        list(
            ecdf = pecdf,
            statistic = actual_stats,
            # we test a greater alternative hypothesis
            # `pecdf(actual_stats)` gives the Pvalue with a permuation stats
            # lower than the actual_stats.
            pvalue_perm = 1L - min(pecdf(actual_stats) + (1L / perm), 1L)
        )
    }
}

#' Rank-Rank Hypergeometric Overlap Map Dots plot
#'
#' @inheritParams rrho_sig_items
#' @param type One of "normal" or "rank", if "normal", the plot value will be
#'   the list1 and list2 used to calculate Rank-Rank Hypergeometric Overlap Map,
#'   Otherwise, the rank of list1 and list2.
#' @inheritDotParams ggplot2::geom_point
#' @return A ggplot Object
#' @importFrom rlang .data
#' @export
rrho_dots <- function(rrho_obj, type = c("normal", "rank"), ...) {
    assert_rrho(rrho_obj)
    type <- match.arg(type)
    list1 <- rrho_obj$rrho_data$list1
    list1 <- data.table::data.table(
        .names = names(list1),
        list1 = unname(list1)
    )
    list2 <- rrho_obj$rrho_data$list2
    list2 <- data.table::data.table(
        .names = names(list2),
        list2 = unname(list2)
    )
    plot_data <- merge(list1, list2, by = ".names")
    if (type == "rank") {
        plot_data[, c("list1", "list2") := lapply(.SD, rank),
            .SDcols = c("list1", "list2")
        ]
    }
    data.table::setDF(plot_data)
    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(.data[["list1"]], .data[["list2"]]) # nolint
    ) +
        ggplot2::geom_point(...)
}

# Since we only use `quadrant_idx_list` and `quadrant_sign` when quadrant isn't
# "all", it's not harmful to assign a non-existent value to them when using
# quadrant "all".
rrho_summary_stats <- function(quadrant, quadrant_idx_list, quadrant_sign, hyper_metric_mat) {
    if (!identical(quadrant, "all")) {
        stats <- vapply(quadrant_idx_list, function(quadrant_idx) {
            max(hyper_metric_mat[quadrant_idx] * quadrant_sign,
                na.rm = TRUE
            )
        }, FUN.VALUE = numeric(1L), USE.NAMES = FALSE)
        sum(stats, na.rm = TRUE)
    } else {
        max(abs(hyper_metric_mat), na.rm = TRUE)
    }
}

rrho_get_direction <- function(x) {
    data.table::fifelse(x <= 0L, "down", "up")
}

assert_rrho <- function(x, arg = rlang::caller_arg(x), call = parent.frame()) {
    assert_class(x, "rrho",
        msg = "{.cls rrho} object returned by {.fn run_rrho} function",
        arg = arg, call = call
    )
}
