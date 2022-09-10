
# Function run GAIA  ------------------------------------------------------

gaia_load_cnv <- function(segmentation_dt, markers_list) {
    message("Loading Copy Number Data")
    stopifnot(inherits(segmentation_dt, "data.frame"))
    segmentation_dt <- data.table::as.data.table(segmentation_dt[1:6])
    data.table::setnames(segmentation_dt, c(
        "Sample", "Chr", "Start", "End", "num_markers", "CN"
    ))
    segmentation_dt[
        , Sample = as.integer(factor(Sample, unique(Sample)))
    ]
    marker_range_list <- lapply(
        names(markers_list), function(marker_range_chr) {
            marker_range_data <- markers_list[[marker_range_chr]]
            GenomicRanges::makeGRangesFromDataFrame(
                df = data.frame(
                    chr = marker_range_chr,
                    start = marker_range_data[1, , drop = TRUE],
                    end = marker_range_data[2, , drop = TRUE]
                ),
                seqnames.field = "chr",
                start.field = "start",
                end.field = "end",
                keep.extra.columns = FALSE
            )
        }
    )
    names(marker_range_list) <- names(markers_list)

    res_temp <- lapply(markers_list, function(chr_probe_matr, nrow) {
        return(
            Matrix::Matrix(
                data = 0, nrow = nrow, ncol = dim(chr_probe_matr)[2],
                dimnames = list(
                    as.character(seq_len(nrow)),
                    as.character(seq_len(dim(chr_probe_matr)[2]))
                ),
                sparse = TRUE
            )
        )
    }, nrow = length(samples_id))


    cnv_list <- split(
        segmentation_dt,
        segmentation_dt[[6]],
        drop = TRUE
    )

    cnv_list <- lapply(cnv_list, function(seg_data, region_data) {
        subject_range <- GenomicRanges::makeGRangesFromDataFrame(
            df = seg_data,
            keep.extra.columns = TRUE,
            seqnames.field = "Chr",
            start.field = "Start",
            end.field = "End"
        )
        subject_range_list_by_chr <- split(
            subject_range,
            GenomeInfoDb::seqnames(subject_range),
            drop = TRUE
        )

        lapply(rlang::set_names(names(res_temp)), function(chr) {
            message(".", appendLF = FALSE)

            sample_probe_matrix <- as.matrix(res_temp[[chr]])

            if (!is.null(subject_range_list_by_chr[[chr]])) {
                hits_res <- GenomicRanges::findOverlaps(
                    marker_range_list[[chr]],
                    subject_range_list_by_chr[[chr]],
                    type = "within"
                )

                if (length(hits_res) > 0) {
                    sample_probe_matrix[cbind(
                        subject_range_list_by_chr[[chr]]$Sample[
                            S4Vectors::subjectHits(hits_res)
                        ],
                        S4Vectors::queryHits(hits_res)
                    )] <- 1
                }
            }

            sample_probe_matrix
        })
    })

    message("\nEnd")
    return(cnv_list)
}

doGAIA <- function(cnv_data, markers_data, aberrations = -1,
                   chromosomes = -1, num_iterations = 10,
                   threshold = 0.25, hom_threshold = 0.12,
                   approximation = FALSE) {
    message("\nPrepare Data")

    stopifnot(inherits(cnv_data, "data.frame"))
    stopifnot(inherits(markers_data, "data.frame"))


    # ..common_chr <- intersect(
    #   cnv_data[[2]], markers_data[[2]]
    # )
    # ..common_chr <- ..common_chr[!is.na(..common_chr)]
    #
    # cnv_data <- cnv_data[ cnv_data[[2]] %in% ..common_chr, ]
    # markers_data <- markers_data[ markers_data[[2]] %in% ..common_chr, ]
    # markers_data <- markers_data[
    #   order(markers_data[[2]], markers_data[[3]]),
    # ]

    # browser()
    markers_obj <- gaia::load_markers(
        marker_matrix = markers_data
    )
    cnv_obj <- gaia_load_cnv(
        segmentation_dt = cnv_data,
        markers_list = markers_obj
    )

    message("\nPerforming Data Preprocessing")
    if (chromosomes == -1) {
        chromosomes <- as.numeric(sort(unique(names(markers_obj))))
        chromosomes <- chromosomes[which(!is.na(chromosomes))]
    } else {
        known_chr <- as.numeric(names(markers_obj))
        if ((length(known_chr[chromosomes]) != length(chromosomes)) ||
            (sum(is.na(known_chr[chromosomes])) > 0)) {
            error_string <- "Error in the list of chromosomes passed as argument.\n"
            error_string <- cat(
                error_string, "The list of the chromosomes that can be analyzed follows:\n",
                known_chr, "\n"
            )
            stop(error_string, call. = FALSE)
        }
    }

    if (aberrations == -1) {
        aberrations <- as.numeric(names(cnv_obj))
        names(aberrations) <- names(cnv_obj)
        if (length(aberrations) > 0 && aberrations[1] == 0) {
            aberrations <- aberrations + 1
        }
    } else {
        aberrations <- as.numeric(names(cnv_obj))
        names(aberrations) <- names(cnv_obj)
        if (length(aberrations) > 0 && aberrations[1] == 0) {
            aberrations <- aberrations + 1
        }
        known_aberr <- as.numeric(names(cnv_obj))
        if ((length(known_aberr[aberrations]) != length(aberrations)) ||
            (sum(is.na(known_aberr[aberrations]) > 0))) {
            error_string <- "Error in the list of aberrations passed as argument.\n"
            error_string <- cat(
                error_string, "The aberrations that can be analyzed follow:\n",
                known_aberr, "\n"
            )
            stop(error_string, call. = FALSE)
        }
    }

    discontinuity <- list()
    if (length(aberrations) == 2 && hom_threshold >= 0) {
        message("\nDone")
        message("\nComputing Discontinuity Matrix")
        for (i in seq_along(chromosomes)) {
            message(".", appendLF = FALSE)
            tmp <- cnv_obj[[2]][[chromosomes[i]]] - cnv_obj[[1]][[chromosomes[i]]]
            tmp_vec <- 0 * c(1:(ncol(tmp) - 1))
            for (k in 1:(ncol(tmp) - 1)) {
                for (z in 1:nrow(tmp)) {
                    tmp_vec[k] <- tmp_vec[k] + abs(tmp[z, k] -
                        tmp[z, k + 1])
                }
            }
            discontinuity[[chromosomes[i]]] <- tmp_vec / nrow(cnv_obj[[2]][[chromosomes[i]]])
        }
    } else {
        if (length(aberrations) != 2 && hom_threshold >= 0) {
            message("\nHomogeneous cannot be applied on the data (data must contain exactly two different kinds of aberrations)\n")
        }
        for (i in seq_along(chromosomes)) {
            tmp <- cnv_obj[[1]][[chromosomes[i]]] - cnv_obj[[1]][[chromosomes[i]]]
            tmp_vec <- 0 * c(1:(ncol(tmp) - 1))
            discontinuity[[chromosomes[i]]] <- tmp_vec
            hom_threshold <- -1
        }
    }
    message("\nDone")


    null_hypothesis_list <- list()
    null_hypothesis_chromosome_list <- list()
    message("Computing Probability Distribution")
    for (k in seq_along(aberrations)) {
        null_hypothesis_chromosome_list <- list()
        for (i in seq_along(chromosomes)) {
            aberrations_index <- (aberrations[k])
            chromosome_index <- chromosomes[i]
            obs_data <- cnv_obj[[aberrations_index]][[chromosome_index]]
            message(".", appendLF = FALSE)
            if (approximation == FALSE) {
                null_hypothesis_chromosome_list[[chromosome_index]] <- gaia::generate_null_hypothesis(
                    obs_data,
                    num_iterations
                )
            }
            if (approximation == TRUE) {
                null_hypothesis_chromosome_list[[chromosome_index]] <- gaia::generate_approx_null_hypothesis(
                    obs_data,
                    num_iterations
                )
            }
        }
        null_hypothesis_list[[aberrations_index]] <- null_hypothesis_chromosome_list
    }
    pvalue_distribution_list <- list()
    for (k in 1:length(aberrations)) {
        tmp_pvalue_list <- list()
        for (i in 1:length(chromosomes)) {
            message(".", appendLF = FALSE)
            tmp_pvalue_list[[chromosomes[i]]] <- rev(cumsum(rev(null_hypothesis_list[[aberrations[k]]][[chromosomes[i]]])))
        }
        pvalue_distribution_list[[aberrations[k]]] <- tmp_pvalue_list
    }
    message("\nDone")


    pvalues_list <- list()
    message("Assessing the Significance of Observed Data")
    for (k in 1:length(aberrations)) {
        pvalue_chromosome_list <- list()
        for (i in 1:length(chromosomes)) {
            message(".", appendLF = FALSE)
            aberrations_index <- (aberrations[k])
            chromosome_index <- chromosomes[i]
            curr_pvalue <- pvalue_distribution_list[[aberrations_index]][[chromosome_index]]
            curr_obs_data <- cnv_obj[[aberrations_index]][[chromosome_index]]
            obs_freq <- apply(curr_obs_data, 2, sum)
            pvalue_chromosome_list[[chromosome_index]] <- round(
                curr_pvalue[obs_freq[] + 1], 7
            )
        }
        pvalues_list[[aberrations_index]] <- pvalue_chromosome_list
    }
    message("\nDone")


    # write file for IGV tools------------------

    # if (length(aberrations) == 2) {
    #   gistic_label <- c("Del", "Amp")
    #   message("Writing ", paste(output_file_name, ".igv.gistic",
    #                             sep = ""), " File for Integrative Genomics Viewer (IGV) Tool")
    #   gistic <- matrix(nrow = 0, ncol = 8)
    #   colnames(gistic) <- c("Type", "Chromosome",
    #                         "Start", "End", "q-value", "G-score",
    #                         "average amplitude", "frequency")
    #   for (k in 1:length(aberrations)) {
    #     tmp_pvalue_list <- list()
    #     for (i in 1:length(chromosomes)) {
    #       message(".", appendLF = FALSE)
    #       curr_qval <- as.numeric(pvalues_list[[aberrations[k]]][[chromosomes[i]]])
    #       curr_qval <- qvalue(curr_qval)
    #       start <- 1
    #       for (z in 2:(length(curr_qval))) {
    #         if (curr_qval[z - 1] != curr_qval[z]) {
    #           end <- z - 1
    #           print_qval <- curr_qval[z - 1]
    #           gistic <- rbind(gistic, c(gistic_label[k],
    #                                     chromosomes[i], markers_obj[[chromosomes[i]]][1,
    #                                                                                   start], markers_obj[[chromosomes[i]]][2,
    #                                                                                                                         end], print_qval, 0, 0, 0))
    #           start <- z
    #         }
    #       }
    #       end <- length(curr_qval)
    #       print_qval <- curr_qval[end]
    #       gistic <- rbind(gistic, c(gistic_label[k], chromosomes[i],
    #                                 markers_obj[[chromosomes[i]]][1, start], markers_obj[[chromosomes[i]]][2,
    #                                                                                                        end], print_qval, 0, 0, 0))
    #     }
    #   }
    #   write.table(gistic, paste(output_file_name, ".igv.gistic",
    #                             sep = ""), sep = "\t", col.names = TRUE,
    #               row.names = FALSE, eol = "\n", quote = FALSE)
    #   message("\nDone")
    # }

    gc()
    if (hom_threshold >= 0) {
        message(
            "Running Homogeneous peel-off Algorithm With Significance Threshold of ",
            threshold, " and Homogenous Threshold of ",
            hom_threshold
        )
    } else {
        message(
            "Running Standard peel-off Algorithm With Significance Threshold of ",
            threshold
        )
    }
    significant_regions_list <- gaia::peel_off(
        pvalues_list, threshold,
        chromosomes, aberrations,
        discontinuity, hom_threshold
    )

    results <- gaia::write_significant_regions(
        markers_obj, significant_regions_list,
        output_file_name = "", chromosomes, aberrations
    )
    return(results)
}
