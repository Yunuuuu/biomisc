test_that("rrho works", {
    set.seed(1L)
    n <- 200
    sample1 <- rnorm(n)
    sample2 <- rnorm(n)
    names(sample1) <- names(sample2) <- paste0("gene", seq_len(n))
    out <- run_rrho(sample1, sample2, 1L)
    testthat::expect_s3_class(out, class = "rrho")
    testthat::expect_snapshot(out)
    testthat::expect_snapshot_output(print(out))

    out2 <- rrho_sig_items(out)
    for (i in seq_along(out2)) {
        testthat::expect_s3_class(out2[[i]], class = "rrho_sig")
    }
    testthat::expect_snapshot(out2)
    testthat::expect_snapshot_output(print(out2))
    testthat::expect_snapshot(rrho_sig_spot(out))
    testthat::expect_snapshot_value(rrho_dots(out), "serialize")
    testthat::expect_snapshot(rrho_correct_pval(out))
})
