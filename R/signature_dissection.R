#' Expectation Maximisation (EM) to divvy up hdp signatures into known cosimic
#' signatures
#'
#' the expectation maximization (EM) algorithm was used to identify pairs of
#' signatures that might explain the observed signature.
#'
#' @details
#' i) For known signatures (e.g. when HDP calls something as “SBS 1”)
#'   - we calculate the cosine similarity between the HDP version of the
#'   signature, and the PCAWG version of the signature.
#'   - If the cosine similarity is below a certain threshold (0.95 seems to work
#'   well), we then apply an expectation maximisation (EM) algorithm to
#'   deconvolute this signature into its constituents. E.g. the preconditioned
#'   HDP SBS1 goes to 38.5% SBS1, 47.5% SBS5, and 15.0% SBS18. Reassuringly, 38%
#'   of the mutations that make up the HDP SBS1 are C to T at CpG, so that fits
#'   well.
#'
#' ii) For novel signatures (i.e. ones that Nicola’s algorithm calls as novel
#' e.g. "SBS N1")
#'   - We always apply EM to break an HDP signature down into a composite of
#'   PCAWG signatures.
#'   - We then reconstitute the signature by adding up the PCAWG signatures
#'   in the proportion that EM gives us (e.g. IDN1 is 27% ID1, 14% ID2, 20% ID5
#'   etc). I have tried only using signatures that contribute >10% of the
#'   mutations to avoid overfitting.
#'   - We compute the cosine similarity of the reconstituted signature to the
#'   known signature. If that is less than a certain threshold, then we consider
#'   that the signature is truly novel. Otherwise, we break it down into its
#'   constituents and present the data in that way.
#'
#' @param signature A numeric vector with the same length of `ncol(targetes)`
#' @param targets A signature matrix with sinature in row and data items in
#'  column (othen COSMIC).
#' @param cos_sim_threthold If the cosine similarity is below a certain
#'   threshold (0.95 seems to work well), we then apply an expectation
#'   maximisation (EM) algorithm to deconvolute this signature into its
#'   constituents.
#' @param emfrac_threshold The minimal fraction returned by EM to rerun EM for
#'  dessection again.
#' @param maxiter An integer specifying the maximal iteration times
#' @param alpha_threshold The minimal error threthold, smaller, the results will
#'  be more accurate.
#' @param unique A scalar logical indicates whether significant (higher than
#'  `cos_sim_threthold`) should be unique, we should set to `TRUE` for priors
#'  signature.
#' @seealso
#' <https://github.com/McGranahanLab/HDP_sigExtraction>
#' @references
#' Frankell, A.M., Dietzen, M., Al Bakir, M. et al. The evolution of lung cancer
#' and impact of subclonal selection in TRACERx. Nature 616, 525–533 (2023).
#' https://doi.org/10.1038/s41586-023-05783-5
#' @export
run_signature_dissection <- function(signature, targets, cos_sim_threthold = 0.9, emfrac_threshold = 0.1, maxiter = 1000L, alpha_threshold = 1e-5, unique = FALSE) {
    assert_class(signature, is.numeric, "numeric")
    assert_class(targets, is.matrix, "matrix")
    if (anyNA(signature)) {
        cli::cli_warn(c(
            "Found {NA} value in {.arg signature}",
            i = "will replace {NA} with 0"
        ))
        signature[is.na(signature)] <- 0
    }
    if (length(signature) != ncol(targets)) {
        cli::cli_abort("{.arg signature} must have the same length of {.code ncol(targets)}")
    }
    if (is.null(rownames(targets))) {
        rownames(targets) <- paste0("Signature", seq_len(nrow(targets)))
    }
    signature_fraction <- run_em(signature, targets,
        maxiter = maxiter, threshold = alpha_threshold
    )
    cos_sim_out <- vapply(seq_len(nrow(targets)), function(i) {
        cos_sim(targets[i, , drop = TRUE], signature)
    }, numeric(1L))
    is_sig <- which(cos_sim_out > cos_sim_threthold)
    nsig <- length(is_sig)
    if ((unique && nsig == 1L) || (!unique && nsig > 0L)) {
        cos_sim_idx <- which.max(cos_sim_out)
        target <- rownames(targets)[cos_sim_idx]
        subset_fraction <- data.table::data.table(
            target = target, fraction = 1L,
            cos_sim = cos_sim_out[cos_sim_idx]
        )
    } else {
        # re-run EM only with signatures with fraction greater than
        # EMfrac_threshold or cosine similarity > 0.8
        constit <- which(signature_fraction > emfrac_threshold)
        constit <- names(constit)
        if (length(constit) <= 1L) {
            constit <- rownames(targets)[
                order(signature_fraction, decreasing = TRUE)
            ][1:2]
        }
        cosine_sigs <- rownames(targets)[cos_sim_out > 0.8] # nolint
        sig_pairs <- utils::combn(unique(c(constit, cosine_sigs)), 2L)
        pairwise_test <- apply(sig_pairs, 2L, function(i) {
            sigs <- targets[i, , drop = FALSE]
            results <- run_em(signature, sigs,
                maxiter = maxiter, threshold = alpha_threshold
            )
            reconstr <- (results[1] * sigs[1L, , drop = TRUE]) +
                (results[2] * sigs[2L, , drop = TRUE])
            data.table::data.table(
                target = names(results),
                fraction = results,
                cos_sim = c(cos_sim(c(reconstr), signature))
            )
        }, simplify = FALSE)
        cos_sim_out <- vapply(pairwise_test, function(x) {
            x$cos_sim[[1L]]
        }, numeric(1L))
        idx <- which.max(cos_sim_out)
        subset_fraction <- pairwise_test[[idx]]
    }
    subset_fraction <- subset_fraction[, !"cos_sim"]
    # reconstitute HDP signatures using the identified COSMIC sigs and calculate
    # the cosine similarity between the original and the reconstituted sig.
    # --> skip if only one COMSIC signature is assigned
    if (length(subset_fraction$target) == 1L && unique) {
        reconsbs <- NULL
        cosine_similarity <- NULL
    } else {
        reconsbs <- .mapply(function(target, fraction) {
            targets[target, , drop = TRUE] * fraction
        }, subset_fraction, NULL)
        reconsbs <- Reduce(`+`, reconsbs)
        cosine_similarity <- cos_sim(reconsbs, signature)
    }
    list(
        signature_fraction = signature_fraction,
        subset_fraction = subset_fraction, reconsbs = reconsbs,
        cosine_similarity = cosine_similarity
    )
}

#' Expectation Maximization (EM) algorithm
#'
#' the expectation maximization (EM) algorithm was used to identify pairs of
#' signatures that might explain the observed signature.
#'
#' @noRd
run_em <- function(signature, targets, maxiter = 1000L, threshold = 1e-5) {
    num_signatures <- nrow(targets)
    # Random start (seems to give ~identical results)
    alpha <- stats::runif(num_signatures)
    alpha <- alpha / sum(alpha)
    # alpha = rep(1/num_signatures,num_signatures) # Uniform start

    for (iter in 1:maxiter) {
        contr <- t(array(alpha, dim = c(num_signatures, ncol(targets)))) * t(targets)
        probs <- contr / array(rowSums(contr), dim = dim(contr))
        probs[is.na(probs)] <- 0
        probs <- probs * signature
        old_alpha <- alpha
        alpha <- colSums(probs) / sum(probs)
        if (sum(abs(alpha - old_alpha)) < threshold) {
            break
        }
    }
    alpha
}
