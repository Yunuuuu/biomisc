#' Estimating CCF with bootstrap
#'
#' @seealso [run_ccf]
#' @note Just for internal usage.
#' @return Modify data in place
#' @noRd
estimate_btstr_ccf <- function(mut_cn_data, sample_field = NULL) {
    assert_pkg("sequenza")
    assert_pkg("boot")
    # check arguments firstly
    assert_df_with_columns(mut_cn_data, c(
        sample_field, "purity", "alt_counts", "ref_counts",
        "minor_cn", "major_cn", "normal_cn"
    ))
    # In order to estimate whether mutations were clonal or subclonal, and the
    # clonal structure of each tumor, a modified version of PyClone was used.
    # For each mutation, two values were calculated, obsCCF and phyloCCF. obsCCF
    # corresponds to the observed cancer cell fraction (CCF) of each mutation.
    # Conversely, phyloCCF corresponds to the phylogenetic CCF of a mutation. To
    # clarify the difference between these two values, consider a mutation
    # present in every cancer cell within a tumor. A subclonal copy number event
    # in one tumor region may lead to loss of this mutation in a subset of
    # cancer cells.  While, the obsCCF of this mutation is therefore below 1,
    # from a phylogenetic perspective the mutation can be considered ‘clonal’ as
    # it occurred on the trunk of the tumor’s phylogenetic tree, and, as such,
    # the phyloCCF may be 1.
    out <- mut_cn_data # modify data in place
    out[, obsVAF := alt_counts / (alt_counts + ref_counts)]

    # estimate the likelihood relating to which copy number the mutation has
    if (is.null(sample_field)) {
        out[, Mt := {
            # for each sample, all mutation share the same types
            types <- new_mut_types(max(minor_cn + major_cn), p = purity[1L])
            get_mt_likelihood(
                obsVAF,
                depths = alt_counts + ref_counts,
                types = types, CNts = minor_cn + major_cn,
                CNns = normal_cn, major_cns = major_cn
            )
        }]
    } else {
        out[, Mt := {
            # for each sample, all mutation share the same types
            types <- new_mut_types(max(minor_cn + major_cn), p = purity[1L])
            get_mt_likelihood(
                obsVAF,
                depths = alt_counts + ref_counts,
                types = types, CNts = minor_cn + major_cn,
                CNns = normal_cn, major_cns = major_cn
            )
        }, by = c(sample_field)]
    }

    # calculate mut_multi
    out[
        ,
        c("mut_multi", "mut_multi_lower", "mut_multi_higher") := calculate_ccf(
            alt_counts, ref_counts,
            CNts = minor_cn + major_cn,
            purity, observed_vafs = obsVAF,
            expected_vafs = NULL, CNns = normal_cn
        )
    ]

    # bootstrap mut_multi CI
    out[
        ,
        c("mut_multi_btstr_lower", "mut_multi_btstr_higher") := bootstrap_cf(
            alt_counts, ref_counts,
            purity = purity, CNts = minor_cn + major_cn, CNns = normal_cn
        )
    ]

    # calculate CCF
    out[,
        paste0("CCF", c("", "_lower", "_higher", "_btstr_lower", "_btstr_higher")) := lapply(.SD, function(x) {
            x / Mt
        }),
        .SDcols = c("mut_multi", "mut_multi_lower", "mut_multi_higher", "mut_multi_btstr_lower", "mut_multi_btstr_higher")
    ]
    out[]
}

utils::globalVariables(c(
    "types_CNn", "types_CNt", "types_Mt", "Mt", "mufreq"
))

# Let's create some functions that can estimate whether early or late
get_mt_likelihood <- function(vafs, depths, types, CNts, CNns = 2L, major_cns) {
    out_list <- .mapply(
        function(vaf, depth, CNt, CNn, major_cn, i) {
            type <- types[types_CNn == CNn & # nolint
                types_CNt == CNt & # nolint
                types_Mt <= major_cn] # nolint
            l <- mufreq_dpois(mufreq = vaf, type$types_mufreq, depth = depth)
            out <- data.table::data.table(l = l / sum(l), Mt = type$types_Mt)
            if (is.na(out$l[1L])) {
                NA_real_
            } else {
                out[which.max(l), Mt] # nolint
            }
        },
        list(
            vaf = vafs, depth = depths, CNt = CNts,
            CNn = rep_len(CNns, length(vafs)), major_cn = major_cns,
            i = seq_along(vafs)
        ), NULL
    )
    unlist(out_list, recursive = FALSE, use.names = FALSE)
}

mufreq_dpois <- function(mufreq, mufreq.model, depth, seq.errors = 0.01, ...) {
    mufreq.model[mufreq.model == 0L] <- seq.errors
    stats::dpois(
        x = round(mufreq * depth),
        lambda = mufreq.model * depth, ...
    )
}

new_mut_types <- function(max_cn, p) {
    assert_class(p, is_scalar_numeric, "scalar {numeric}")
    types <- lapply(1:2, function(CNn, max_cn) {
        sequenza::mufreq.types.matrix(CNt.min = 1L, CNt.max = max_cn, CNn = CNn)
    }, max_cn = max_cn)
    types <- data.table::rbindlist(types)
    types <- types[Mt >= 1L] # nolint
    # theoretical.mufreq returns the theoretical mutation frequency at a single
    # specific position, given values of cellularity, copy number in the normal
    # and tumor samples at that position, and the number of mutated alleles.
    types[, mufreq := mapply(function(CNt, CNn, Mt, p) { # nolint
        sequenza::theoretical.mufreq(
            CNt = CNt, Mt = Mt, cellularity = p, CNn = CNn # nolint
        )
    }, CNn = CNn, CNt = CNt, Mt = Mt, MoreArgs = list(p = p))]
    data.table::setnames(types, function(x) paste0("types_", x))
    types
}

# should run for each mutation
bootstrap_cf <- function(alt_counts, ref_counts, purity, CNts, CNns, times = 1000L) {
    out_list <- .mapply(function(alt_count, ref_count, p, CNn, CNt, times) {
        if (ref_count == 0L) {
            ci_list <- calculate_ccf(
                alt_counts = alt_count,
                ref_counts = ref_count,
                CNts = CNt, p,
                observed_vafs = alt_count / (alt_count + ref_count),
                expected_vafs = NULL
            )[2:3]
            unlist(ci_list, use.names = FALSE)
        } else {
            x <- c(rep(1L, alt_count), rep(0L, ref_count))
            theta <- function(x, i) {
                data <- x[i]
                calculate_mut_multi(CNt, sum(data) / length(data), p, CNn)
            }
            boot_res <- boot::boot.ci(
                boot::boot(x, theta, R = times),
                type = "norm"
            )
            c(boot_res$normal[2L], boot_res$normal[3L])
        }
    }, list(
        alt_count = alt_counts, ref_count = ref_counts,
        p = purity, CNt = CNts, CNn = CNns
    ), list(times = times))
    data.table::transpose(out_list)
}
