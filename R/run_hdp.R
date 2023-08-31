#' Posterior sampling chain across activated DPs.
#'
#' Run a Gibbs sampler over the activated DP nodes of a Hierarchichal Dirichlet
#' Process.  Each iteration re-assigns the cluster allocation of every data
#' item.  Run `burnin` iterations, and then collect `n` samples from the chain
#' with `space` iterations between each collected sample.
#'
#' @param matrix A matrix of counts with one row for every sample (same order
#'   as dpindex) and one column for every data category.
#' @param priors A matrix of prior distributions (columns must each sum to 1,
#'   number of rows matches number of data categories). If not NULL, the number
#'   of rows must match the number of columns in matrix.
#' @param prior_pseudoc Vector of pseudocounts contributed by each prior
#'   distribution. Only used when priors are provided. The frozen nodes
#'   signifying prior components contain some number of pseudo-data items - use
#'   these pseudocounts to specify the relative weighting of each prior
#'   component compared to the volume of observed data. If NULL, will be
#'   `rep_len(1000L, ncol(priors))`.
#' @param dp_tree A data.frame specifying the children-parent relationships. The
#'   first column must match `rownames(matrix)`, and the second define the
#'   parent, the third difine the grandparent and more ancestors.
#' @param initcc Number of data clusters to start with (every data item is
#'   randomly assigned to a cluster to start with). Will coerced to integer. See
#'   [dp_activate][hdp::dp_activate].
#' @param n_posterior Specifying the number of independent posterior sampling
#'   chains to run. Will coerced to integer.
#' @inheritDotParams hdp::hdp_posterior -hdp -seed
#' @param seed The seed that can be set to reproduce output.
#' @seealso [hdp_posterior][hdp::hdp_posterior]
#' @return A `HDP` object.
#' @examples
#' run_hdp(
#'     matrix = hdp::example_data_hdp_prior,
#'     priors = hdp::example_known_priors,
#'     burnin = 10L, # 10000L
#'     n = 10L, # 100L
#'     space = 200L
#' )
#' @export
run_hdp <- function(
    matrix, priors = NULL, prior_pseudoc = NULL,
    dp_tree = NULL, initcc = 10L, ..., n_posterior = 15L, seed = 1234L) {
    assert_pkg("hdp")
    assert_class(matrix, is.matrix, "matrix")
    if (anyNA(matrix)) {
        cli::cli_abort("NA is not allowed in {.arg matrix}")
    }
    assert_class(priors, is.matrix, "matrix", null_ok = TRUE)
    if (!is.null(priors)) {
        if (is.null(colnames(matrix))) {
            if (ncol(matrix) != nrow(priors)) {
                cli::cli_abort("{.code ncol(matrix)} and {.code nrow(priors)} must be equal")
            }
        } else {
            if (!all(colnames(matrix) == rownames(priors))) {
                cli::cli_abort("{.code colnames(matrix)} and {.code rownames(priors)} must be the same")
            }
        }
    }
    assert_class(initcc, is_scalar_numeric, "scalar numeric", cross_msg = NULL)
    assert_class(n_posterior, is_scalar_numeric, "scalar numeric",
        cross_msg = NULL
    )
    assert_class(seed, is_scalar_numeric, "scalar numeric", cross_msg = NULL)
    oldseed <- get0(".Random.seed", envir = .GlobalEnv)
    if (is.null(oldseed)) {
        on.exit(rm(".Random.seed", envir = .GlobalEnv))
    } else {
        on.exit(assign(".Random.seed", value = oldseed, envir = .GlobalEnv))
    }
    set.seed(seed,
        kind = "Mersenne-Twister",
        normal.kind = "Inversion",
        sample.kind = "Rejection"
    )
    initcc <- as.integer(initcc)
    n_posterior <- as.integer(n_posterior)
    if (is.null(rownames(matrix))) {
        rownames(matrix) <- seq_len(nrow(matrix))
    }
    dp_tree <- hdp_prepare_tree(dp_tree, matrix)

    ##################################
    ###  initialise HDP structure  ###
    ##################################
    cli::cli_inform("Preparing data")
    if (!is.null(priors)) {
        ### with priors ###
        nps <- ncol(priors)
        assert_class(prior_pseudoc, is.numeric, "numeric",
            null_ok = TRUE, cross_msg = NULL
        )
        assert_length(prior_pseudoc, nps, scalar_ok = TRUE, null_ok = TRUE)
        # (don’t activate the frozen pseudo-count nodes for the prior
        # signatures)
        initial <- nps + 2L
        initcc <- initcc + nps
        hdp_data <- hdp::hdp_prior_init(
            prior_distn = priors,
            prior_pseudoc = prior_pseudoc %||% rep_len(1000L, nps),
            hh = rep_len(1, ncol(matrix)), # must not a integer
            # all priors share the same hyperparameters
            alphaa = rep_len(1L, 2L),
            alphab = rep_len(1L, 2L)
        )

        ## add a parent none-prior node
        hdp_data <- hdp::hdp_addconparam(hdp_data,
            alphaa = 1L, alphab = 1L
        )
        hdp_data <- hdp::hdp_adddp(hdp_data,
            numdp = 1L, ppindex = 1L,
            cpindex = 1L + max(hdp::cpindex(hdp_data))
        )
    } else {
        initial <- 1L
        hdp_data <- hdp::hdp_init(
            ppindex = 0L, cpindex = 1L,
            hh = rep_len(1, ncol(matrix)), # must not a integer
            alphaa = 1L, alphab = 1L
        )
    }

    # add dp and conparam --------------------------
    number_groups <- ncol(dp_tree)
    for (i in rev(seq_len(number_groups))) {
        numdp <- length(unique(dp_tree[[i]]))
        if (i == number_groups) {
            idx <- rep_len(1L, numdp)
            # the first parent idx should run over the length of current ppindex
            # like in some ways there are some priors, but shouldn't run over
            # the unique number
            ppidx <- rep_len(length(hdp::ppindex(hdp_data)), numdp)
        } else {
            # only include the first parent in the next column
            # Others are ignored
            idx <- dp_tree[
                !duplicated(dp_tree[[i]]),
                factor(.SD[[1L]], unique(.SD[[1L]])),
                .SDcols = i + 1L
            ]
            idx <- as.integer(idx)
            ppidx <- idx + max(hdp::ppindex(hdp_data))
        }
        # all children with the same parent should have same concentration
        num_conparam <- length(unique(idx))
        hdp_data <- hdp::hdp_addconparam(hdp_data,
            alphaa = rep_len(1L, num_conparam),
            alphab = rep_len(1L, num_conparam)
        )
        hdp_data <- hdp::hdp_adddp(hdp_data,
            numdp = numdp,
            ppindex = ppidx,
            cpindex = idx + max(hdp::cpindex(hdp_data))
        )
    }

    # add data to leaf nodes (one per sample, in row order of matrix)
    terminal <- hdp::numdp(hdp_data)
    hdp_data <- hdp::hdp_setdata(hdp_data,
        dpindex = (terminal - nrow(matrix) + 1L):terminal,
        matrix
    )

    ############
    ### run posterior sampling chains
    ############
    cli::cli_inform("Running chain sampling")
    posteriors <- lapply(seq_len(n_posterior), function(i) {
        # with different dp_activate seeds to start the clustering from a
        # different random starting point each time.
        new_hdp_data <- hdp::dp_activate(
            hdp = hdp_data,
            dpindex = initial:terminal,
            initcc = initcc
        )
        hdp::hdp_posterior(hdp = new_hdp_data, ...)
    })

    structure(
        list(
            input = matrix, priors = priors,
            posteriors = posteriors,
            components = NULL, statistics = NULL
        ),
        class = "HDP"
    )
}

hdp_prepare_tree <- function(dp_tree, matrix, arg1 = rlang::caller_arg(dp_tree), arg2 = rlang::caller_arg(matrix), call = parent.frame()) {
    assert_class(dp_tree, "data.frame", null_ok = TRUE, arg = arg1, call = call)
    if (!is.null(dp_tree)) {
        dp_tree <- data.table::as.data.table(dp_tree)
        if (!all(dp_tree[[1L]] == rownames(matrix))) {
            cli::cli_abort(
                "The first column of {.arg {arg1}} must match {.code rownames({arg2})}",
                call = call
            )
        }
        if (ncol(dp_tree) >= 2L) {
            all_items <- names(dp_tree)
            for (kk in seq_len(ncol(dp_tree) - 1L)) {
                assert_nest(dp_tree, all_items[[kk + 1L]], all_items[[kk]])
            }
        }
    } else {
        dp_tree <- data.table::data.table(sample = rownames(matrix))
    }
    dp_tree
}

#' @export
print.HDP <- function(x, ...) {
    cli::cli_text("A {.cls HDP} object with {.val {length(x$posteriors)}} Posterior sampling chain{?s}")
    if (!is.null(x$statistics)) {
        stats <- x$statistics
        comps <- sort(rownames(stats$signif_signatures))
        new_components <- comps[startsWith(comps, "N")]
        prior_components <- comps[startsWith(comps, "P")]
        components_to_nsamples <- colSums(stats$signif_exposures > 0L)
        components_to_counts <- colSums(stats$signif_exposures_counts)

        if (length(prior_components)) {
            values <- paste(
                "A total of {.val {",
                components_to_counts[prior_components],
                "}} count{?s} in {.val {",
                components_to_nsamples[prior_components],
                "}} sample{?s}"
            )
            names(values) <- prior_components
            cli::cli_text("Components in Priors")
            cli::cli_dl(values)
        }
        if (length(new_components)) {
            values <- paste(
                "A total of {.val {",
                components_to_counts[new_components],
                "}} count{?s} in {.val {",
                components_to_nsamples[new_components],
                "}} sample{?s}"
            )
            names(values) <- new_components
            cli::cli_text("New identified Components")
            cli::cli_dl(values)
        }
    }
    invisible(x)
}

#' @export
`$.HDP` <- function(x, i) {
    NextMethod()
}

#' @export
`$<-.HDP` <- function(x, i, value) {
    NextMethod()
}

#' @export
`[[.HDP` <- function(x, i) {
    NextMethod()
}

#' @export
`[[<-.HDP` <- function(x, i, value) {
    NextMethod()
}

#' @export
`[.HDP` <- function(x, i) {
    NextMethod()
}

#' @export
`[<-.HDP` <- function(x, i, value) {
    NextMethod()
}

#############################################################################
#' Extract hdp results
#'
#' This involves collecting multiple independent HDP sampling chains into an
#' hdpSampleMulti object using [hdp_multi_chain][hdp::hdp_multi_chain].
#' Components are then extracted using
#' [hdp_extract_components][hdp::hdp_extract_components].
#' @param x A `HDP` object returned by [run_hdp].
#' @param ...
#'  * dpindices: Indices of DP nodes to extract. Default: `NULL`.
#'  * sig_active_cutoff: A numeric of the minimal weight to regard a component
#'  as active. Default: `0.1`.
#'  * remove_zero_lower_ci: A scalar logical indicates exposures with
#'    zero lower confidence interval should be removed.
#'  * cohort_threshold: A numeric of the minimal proportion (if <1L) or number
#'    of samples (if >= 1L) to regard a component as active. Default: `0.05`.
#' @return A `HDP` object with `components` and `statistics` added.
#' @export
hdp_data <- function(x, ...) {
    assert_pkg("hdp")
    assert_class(x, "HDP", "{.cls HDP} object returned by {.fn run_hdp}",
        cross_msg = NULL
    )
    hdp_multi_chain <- hdp::hdp_multi_chain(x$posteriors)
    x$components <- hdp::hdp_extract_components(hdp_multi_chain)
    x$statistics <- hdp_data_internal(
        hdpsample = x$components,
        input_matrix = x$input, ...
    )
    x
}

hdp_data_internal <- function(
    hdpsample, input_matrix, dpindices = NULL,
    sig_active_cutoff = 0.1, remove_zero_lower_ci = TRUE,
    cohort_threshold = 0.05) {
    assert_class(sig_active_cutoff, function(x) {
        is_scalar_numeric(x) && data.table::between(x, 0L, 1L)
    }, "scalar {.cls numeric} in [0, 1]", NULL)
    assert_class(
        remove_zero_lower_ci, rlang::is_scalar_logical,
        "scalar {.cls logical}", NULL
    )
    assert_class(cohort_threshold, function(x) {
        is_scalar_numeric(x) && x >= 0L
    }, "scalar {.cls numeric} not less than 0", NULL)
    if (length(hdp::comp_categ_counts(hdpsample)) == 0L) {
        cli::cli_abort("No component info for hdpsample. First run hdp_extract_components")
    }
    # define signature
    comp_distn <- hdp::comp_categ_distn(hdpsample)
    signatures <- comp_distn$mean
    if (ncol(input_matrix) != ncol(signatures)) {
        cli::cli_abort("{.arg input_matrix} is not compatible with {.arg hdpsample}")
    }
    colnames(signatures) <- colnames(input_matrix)

    # define exposure
    dp_distn <- hdp::comp_dp_distn(hdpsample)
    ndp <- nrow(dp_distn$mean)
    ncomp <- ncol(dp_distn$mean) # nolint
    if (methods::is(hdpsample, "hdpSampleChain")) {
        dps <- hdp::dp(hdp::final_hdpState(hdpsample))
        pps <- hdp::ppindex(hdp::final_hdpState(hdpsample))
    } else if (methods::is(hdpsample, "hdpSampleMulti")) {
        dps <- hdp::dp(hdp::final_hdpState(hdp::chains(hdpsample)[[1]]))
        pps <- hdp::ppindex(hdp::final_hdpState(hdp::chains(hdpsample)[[1]]))
    }
    if (is.null(dpindices)) {
        dpindices <- (length(pps) - nrow(input_matrix) + 1L):length(pps)
    } else if (!is.numeric(dpindices) ||
        !any(round(dpindices) == dpindices) ||
        any(dpindices < 1L) || any(dpindices > ndp)) {
        cli::cli_abort("dpindices must be integers between 1 and {ndp}")
    }
    dps <- dps[dpindices]
    pps <- pps[dpindices]

    # caculate signatures and exposures --------------------------------
    exposures <- dp_distn$mean[dpindices, , drop = FALSE]
    if (nrow(exposures) != nrow(input_matrix)) {
        cli::cli_warn("{.arg dpindices} is not compatible with {.arg input_matrix}")
    } else {
        rownames(exposures) <- rownames(input_matrix)
    }
    reconstructed <- exposures %*% signatures

    # reconstruction error --------------------------------------------
    RMSE <- vapply(seq_len(nrow(input_matrix)), function(i) {
        rmse(input_matrix[i, , drop = TRUE], reconstructed[i, , drop = TRUE])
    }, numeric(1L))
    nRMSE <- RMSE / rowMeans(input_matrix)
    cosineSimilarity <- vapply(seq_len(nrow(input_matrix)), function(i) {
        cos_sim(input_matrix[i, , drop = TRUE], reconstructed[i, , drop = TRUE])
    }, numeric(1L))

    # significant signature -------------------------------------------
    signif_exposures <- exposures
    cis <- dp_distn$cred.int[dpindices]
    if (remove_zero_lower_ci) {
        nonsig <- lapply(cis, function(x) which(x[1L, ] == 0L))
        for (i in seq_along(nonsig)) {
            signif_exposures[i, nonsig[[i]]] <- 0L
        }
    }
    signif_exposures <- signif_exposures[
        , colnames(signif_exposures) != "0",
        drop = FALSE
    ]
    signif_signatures <- signatures[rownames(signatures) != "0", , drop = FALSE]
    signif_metrics <- data.table::as.data.table(signif_exposures,
        keep.rownames = "sample_id"
    )
    signif_metrics <- data.table::melt(signif_metrics,
        id.vars = "sample_id",
        variable.name = "components",
        variable.factor = FALSE,
        value.name = "exposures"
    )
    signif_metrics$sig_active <- signif_metrics$exposures >= sig_active_cutoff
    if (cohort_threshold < 1L) {
        cohort_threshold <- cohort_threshold * nrow(exposures)
    }
    cohort_threshold <- as.integer(cohort_threshold)
    cli::cli_inform(
        "Filtering components in less than {cohort_threshold} sample{?s}"
    )
    signif_metrics[, active_samples := sum(sig_active, na.rm = TRUE), # nolint
        by = "components"
    ]
    signif_metrics[, sig_cohort := active_samples >= cohort_threshold]
    excluded_components <- signif_metrics[(!sig_cohort), unique(components)]
    signif_signatures <- signif_signatures[
        !rownames(signif_signatures) %in% excluded_components, ,
        drop = FALSE
    ]
    signif_exposures <- signif_exposures[
        , !colnames(signif_exposures) %in% excluded_components,
        drop = FALSE
    ]
    mutBurden <- data.table::data.table(
        sample = rownames(input_matrix),
        nMut = rowSums(input_matrix)
    )
    list(
        reconstructed = reconstructed,
        signatures = signatures,
        exposures = exposures,
        numdata = vapply(dps, hdp::numdata, integer(1L)),
        exposures_counts = round(exposures * mutBurden$nMut),
        signif_signatures = signif_signatures,
        signif_exposures = signif_exposures,
        signif_exposures_counts = round(signif_exposures * mutBurden$nMut),
        # signif_metrics = signif_metrics,
        excluded_components = excluded_components,
        RMSE = RMSE, nRMSE = nRMSE,
        cosineSimilarity = cosineSimilarity,
        mutBurden = mutBurden
    )
}

utils::globalVariables(c(
    "active_samples", "sig_active", "sig_cohort", "components"
))

cos_sim <- function(x, y) {
    crossprod(x, y) / sqrt(crossprod(x) * crossprod(y))
}

rmse <- function(actual, predicted) {
    sqrt(mse(actual, predicted))
}

mse <- function(actual, predicted) {
    mean(se(actual, predicted))
}

se <- function(actual, predicted) {
    (actual - predicted)^2L
}
