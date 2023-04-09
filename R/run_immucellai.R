#' Immune Cell Abundance Identifier
#'
#' ImmuCellAI (Immune Cell Abundance Identifier,
#' <http://bioinfo.life.hust.edu.cn/ImmuCellAI/>) is a tool to estimate the
#' abundance of 24 immune cells from gene expression dataset including RNA-Seq
#' and microarray data
#'
#' @param sample_exp Gene expression matrix
#' @param data_type One of "rnaseq" and "microarray"
#' @return a data.frame
#' @examples
#' sample_exp <- readRDS(system.file(
#'     "extdata", "immucellai", "sample_exp.rds",
#'     package = "biomisc"
#' ))
#' run_immucellai(sample_exp, "microarray")
#' @seealso <https://github.com/lydiaMyr/ImmuCellAI>
#' @export
run_immucellai <- function(sample_exp, data_type = c("microarray", "rnaseq")) {
    assert_pkg("GSVA")
    assert_pkg("pracma")
    assert_pkg("quadprog")
    assert_class(
        sample_exp,
        function(x) is.matrix(x) && is.numeric(x),
        msg = "numeric {.cls matrix} object"
    )
    data_type <- match.arg(data_type)

    paper_marker <- readRDS(
        system.file("extdata", "immucellai", "paper_marker.rds",
            package = "biomisc"
        )
    )
    marker_exp <- readRDS(
        system.file("extdata", "immucellai", "marker_exp.rds",
            package = "biomisc"
        )
    )
    # save the expression value of common genes
    common_genes <- intersect(unlist(paper_marker), rownames(marker_exp))
    common_genes <- intersect(common_genes, rownames(sample_exp))

    sam_exp <- sample_exp[common_genes, , drop = FALSE]
    marker_exp <- marker_exp[common_genes, , drop = FALSE]
    if (identical(data_type, "rnaseq")) sam_exp <- log2(sam_exp + 1L)

    # prepare marker tag matrix
    marker_tag_mat <- lapply(paper_marker, function(markers) {
        data.table::fifelse(common_genes %chin% markers, 1L, 0L)
    })

    # prepare sample expression matrix
    sam_exp <- apply(sam_exp, 2L, function(exp) {
        res <- lapply(names(paper_marker), function(cell) {
            exp / marker_exp[, cell, drop = TRUE] * marker_tag_mat[[cell]]
        })
        data.table::setDT(res)
        rowSums(res)
    })
    rownames(sam_exp) <- common_genes

    # run GSVA
    result <- GSVA::gsva(
        sam_exp, paper_marker,
        method = "ssgsea",
        ssgsea.norm = TRUE
    )

    if (ncol(result) < 3L) {
        result[result < 0L] <- 0L
    } else {
        result <- result - apply(result, 1L, min)
    }
    compensation_matrix <- readRDS(
        system.file("extdata", "immucellai", "compensation_matrix.rds",
            package = "biomisc"
        )
    )
    result_norm <- compensation(result, compensation_matrix)
    infiltrating_score <- colSums(
        result_norm[
            c("Bcell", "CD4_T", "CD8_T", "DC", "Macrophage", "Monocyte", "Neutrophil", "NK"), ,
            drop = FALSE
        ]
    )
    infiltrating_score <- (infiltrating_score / max(infiltrating_score)) * 0.9
    result_mat <- rbind(result_norm, infiltrating_score = infiltrating_score)
    res <- data.table::as.data.table(t(result_mat), keep.rownames = "Samples")
    data.table::setDF(res)
    res
}

compensation <- function(raw_score, compensation_matrix) {
    diag(compensation_matrix) <- 1L
    common_cells <- rownames(raw_score)[
        rownames(raw_score) %chin% rownames(compensation_matrix)
    ]
    scores <- apply(raw_score[common_cells, , drop = FALSE], 2L, function(x) {
        pracma::lsqlincon(
            compensation_matrix[common_cells, common_cells], x,
            lb = 0L
        )
    })
    scores[scores < 0L] <- 0L
    rownames(scores) <- common_cells
    scores
}
