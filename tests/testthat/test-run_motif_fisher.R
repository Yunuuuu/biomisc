testthat::test_that("run_motif_fisher works", {
    maf <- readRDS(system.file("extdata", "absolute",
        "run_absolute_example_maf.rds",
        package = "biomisc"
    ))
    ref_genome <- suppressWarnings(BSgenome::getBSgenome("hg19"))
    testthat::expect_snapshot(
        suppressWarnings(run_motif_fisher(
            maf[, c("Tumor_Sample_Barcode", "Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele2")],
            ref_genome = ref_genome
        ))
    )
})
