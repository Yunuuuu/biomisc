test_that("run_signature_dissection works", {
    hdp_test <- readRDS("../testdata/EM_split_signature/test_data.rds")
    hsbs <- hdp_test$hsbs
    psbs <- hdp_test$psbs
    signatures <- t(psbs)
    check_sigs <- colnames(hsbs)
    maxiter_EM <- 1000L
    EMfrac_threshold <- 0.1
    mutations <- hsbs[, check_sigs, drop = FALSE]

    # caculate cosineSimilarity
    cosmat <- data.frame(matrix(nrow = ncol(hsbs), ncol = ncol(psbs)))
    rownames(cosmat) <- colnames(hsbs)
    colnames(cosmat) <- colnames(psbs)

    for (i in seq_len(nrow(cosmat))) {
        for (j in seq_len(ncol(cosmat))) {
            cosmat[i, j] <- cos_sim(
                x = hsbs[, rownames(cosmat)[i]],
                y = psbs[, colnames(cosmat)[j]]
            )
        }
    }
    colnames(cosmat) <- paste0("cosmic_", colnames(cosmat))
    rownames(cosmat) <- paste0("hdp_", rownames(cosmat))
    cosineSim_threshold <- 0.9

    EM_algorithm <- function(output_signature, input_signatures, maxiter = 1000) {
        num_signatures <- nrow(input_signatures)
        alpha <- runif(num_signatures)
        alpha <- alpha / sum(alpha) # Random start (seems to give ~identical results)
        # alpha = rep(1/num_signatures,num_signatures) # Uniform start

        for (iter in 1:maxiter) {
            contr <- t(array(alpha, dim = c(num_signatures, ncol(input_signatures)))) * t(input_signatures)
            probs <- contr / array(rowSums(contr), dim = dim(contr))
            probs[is.na(probs)] <- 0
            probs <- probs * output_signature
            old_alpha <- alpha
            alpha <- colSums(probs) / sum(probs)
            if (sum(abs(alpha - old_alpha)) < 1e-5) {
                break
            }
        }
        return(alpha)
    }

    # test
    for (check_sig in check_sigs) {
        # check_sig <- "P2(SBS5)"
        set.seed(1L)
        signature_fraction2 <- run_em(mutations[, check_sig], signatures)

        set.seed(1L)
        signature_fraction <- lapply(check_sig, function(x) {
            output_signature <- mutations[, x]
            # output_signature[is.na(output_signature)] <- 0
            as.matrix(EM_algorithm(output_signature, signatures, maxiter_EM))
        })
        signature_fraction <- Reduce(cbind, signature_fraction)
        colnames(signature_fraction) <- check_sig
        testthat::expect_equal(signature_fraction[, 1L], signature_fraction2)


        # re-run EM only with signatures with fraction greater than EMfrac_threshold and cosine similarity > 0.8
        subset_signature_fractions <- lapply(check_sig, function(x) {
            print(x)
            # x <- check_sig
            constit <- names(which(signature_fraction[, x] > EMfrac_threshold))
            if (length(sub("cosmic_", "", colnames(cosmat)[which(cosmat[paste0("hdp_", x), ] >= cosineSim_threshold)])) > 0) {
                if (length(sub("cosmic_", "", colnames(cosmat)[which(cosmat[paste0("hdp_", x), ] >= cosineSim_threshold)])) > 1) {
                    cosmic <- cosmat[paste0("hdp_", x), which(cosmat[paste0("hdp_", x), ] >= cosineSim_threshold)]
                    cosmic <- sub("cosmic_", "", colnames(cosmic)[which.max(cosmic)])
                } else {
                    cosmic <- sub("cosmic_", "", colnames(cosmat)[which(cosmat[paste0("hdp_", x), ] >= cosineSim_threshold)])
                }
                out <- data.frame(hdp = x, cosmic = cosmic, fraction = 1)
                return(out)
            }
            cosine_sigs <- sub("cosmic_", "", colnames(cosmat)[which(cosmat[paste0("hdp_", x), ] > 0.8)])

            if (length(constit) <= 1) {
                constit <- names(sort(signature_fraction[, x], decreasing = T)[1:2])
            }

            sig_pairs <- t(combn(unique(c(constit, cosine_sigs)), 2))
            pairwise_test <- lapply(1:nrow(sig_pairs), function(i) {
                sigs <- as.matrix(psbs[, sig_pairs[i, ], drop = F])
                input_signatures <- t(sigs)
                output_signature <- hsbs[, x]
                results <- EM_algorithm(output_signature, input_signatures, maxiter_EM)
                # while(any(results <= EMfrac_threshold)){
                #   constit     <- names(which(results > EMfrac_threshold))
                #   sigs        <- as.matrix(psbs[,constit, drop = F])
                #   input_signatures <- t(sigs)
                #   output_signature <- hsbs[,x]
                #   results    <- EM_algorithm(output_signature, input_signatures, maxiter_EM)
                # }

                reconstr <- results[1] * psbs[, sig_pairs[i, 1]] + results[2] * psbs[, sig_pairs[i, 2]]
                cosSim <- lsa::cosine(reconstr, output_signature)
                out <- data.frame(hdp = x, cosmic = names(results), fraction = as.numeric(results), cosineSimilarity = cosSim)
                return(out)
            })

            cosineSims <- sapply(pairwise_test, function(y) y$cosineSimilarity[1])
            index <- which(cosineSims == max(cosineSims))

            out <- pairwise_test[[index]]
            out <- out[, -1 * ncol(out)]
            return(out)
        })

        set.seed(1L)
        testthat::expect_equal(
            as.data.frame(
                run_signature_dissection(
                    mutations[, check_sig], signatures,
                    unique = FALSE
                )$subset_fraction[, list(
                    cosmic = target, fraction = fraction
                )]
            ),
            subset_signature_fractions[[1L]][c("cosmic", "fraction")],
            tolerance = sqrt(.Machine$double.eps)
        )
    }
})
