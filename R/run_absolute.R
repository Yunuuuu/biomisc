#' Automate ABSOLUTE calling for multiple samples in parallel ways
#'
#' @description This function is based on package
#'   \href{https://github.com/ShixiangWang/DoAbsolute}{DoAbsolute} by adding
#'   parallel evaluation and adjusting personal convention
#'
#' @details \href{https://www.nature.com/articles/nbt.2203}{ABSOLUTE} is a
#'   famous software developed by Broad Institute, Using
#'   install.packages("https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/absolute/ABSOLUTE_1.0.6.tar.gz",
#'   repos = NULL, type = "source") to install it. However, the
#'   \code{\link[ABSOLUTE]{RunAbsolute}} function points to estimate one sample
#'   each time and sets no default values. \code{\link{run_absolute}} helps
#'   users set default parameters based on
#'   \href{https://www.genepattern.org/modules/docs/ABSOLUTE}{ABSOLUTE
#'   documentation} and provides an uniform interface to input data easily and
#'   run \code{\link[ABSOLUTE]{RunAbsolute}} parallelly.
#'
#'   If calling for a sample failed, the error message will be written to
#'   \code{error.log} under subdirectory \code{RunAbsolute} of
#'   \code{results_dir} directory.
#'
#'   More detail about how to analyze ABSOLUTE results please see
#'   \href{https://www.genepattern.org/analyzing-absolute-data}{analyzing-absolute-data}.
#'
#' @param seg a \code{data.frame} containing columns "Chromosome", "Start",
#'   "End", "Num_Probes", "Segment_Mean". If providing multiple samples, seg
#'   should contain a column "Sample" to identify different samples
#' @param maf MAF, default is \code{NULL}, can provided as \code{data.frame}.
#' @param sigma_p Provisional value of excess sample level variance used for
#'   mode search. Default: \code{0}
#' @param max_sigma_h Maximum value of excess sample level variance (Eq. 6).
#'   Default: \code{0.015}
#' @param min_ploidy Minimum ploidy value to consider. Solutions implying lower
#'   ploidy values will be discarded. Default: \code{0.95}
#' @param max_ploidy Maximum ploidy value to consider. Solutions implying
#'   greater ploidy values will be discarded. Default: \code{10}
#' @param primary_disease Primary disease of the sample. Default: \code{NA}
#' @param platform one of \code{"SNP_6.0"}, \code{"Illumina_WES"},
#'   \code{"SNP_250K_STY"}. Default: \code{"SNP_6.0"}
#' @param results_dir directory path used to store result files. Default:
#'   `"ABSOLUTE"`
#' @param max_as_seg_count Maximum number of allelic segments. Samples with a
#'   higher segment count will be flagged as 'failed'. Default: \code{1500}
#' @param max_non_clonal Maximum genome fraction that may be modeled as
#'   non-clonal (subclonal SCNA). Solutions implying greater values will be
#'   discarded. Default: \code{0.05}
#' @param max_neg_genome Maximum genome fraction that may be modeled as
#'   non-clonal with copy-ratio below that of clonal homozygous deletion.
#'   Solutions implying greater values will be discarded. Default: \code{0.005}
#' @param copy_num_type The type of copy number to be handled. Either total or
#'   allelic. total is what this code for. Default: \code{"total"}
#' @param min_mut_af Minimum mutation allelic fraction. Mutations with lower
#'   allelic fractions will be filtered out before analysis. Default: \code{0.1}
#' @author Yun \email{yunyunpp96@@outlook.com}
#' @return  Side effect. \cr \cr All ABSOLUTE called results (see
#'   \code{\link[ABSOLUTE]{RunAbsolute}}) were kept in directory
#'   \code{file.path(results_dir, "RunAbsolute")}. \cr \cr All summarized
#'   results from multiple ABSOLUTE calling (see
#'   \code{\link[ABSOLUTE]{CreateReviewObject}}) were kept in
#'   \code{file.path(results_dir, "CreateReviewObject")}. \cr \cr All reviewed
#'   results (see \code{\link[ABSOLUTE]{ExtractReviewedResults}}) were kept in
#'   \code{ file.path(results_dir, "reviewed")}
#' @examples
#' \donttest{
#' seg <- readRDS(system.file("extdata", "absolute",
#'     "run_absolute_example_seg.rds", package = "biomisc" ))
#' maf <- readRDS(system.file("extdata", "absolute",
#'     "run_absolute_example_maf.rds", package = "biomisc" ))
#' run_absolute(
#'     seg = seg, maf = maf,
#'     results_dir = file.path(tempdir(), "ABSOLUTE")
#' )
#' }
#' @references Carter, S., Cibulskis, K., Helman, E. et al. Absolute
#'   quantification of somatic DNA alterations in human cancer. Nat Biotechnol
#'   30, 413–421 (2012). \url{https://doi.org/10.1038/nbt.2203}
#' @export
run_absolute <- function(seg, maf = NULL, sigma_p = 0, max_sigma_h = 0.015,
                         min_ploidy = 0.95, max_ploidy = 10,
                         primary_disease = NA,
                         platform = NULL,
                         results_dir = "ABSOLUTE",
                         max_as_seg_count = 1500, 
                         max_neg_genome = 0.005,
                         max_non_clonal = 0.05,
                         copy_num_type = NULL,
                         min_mut_af = 0.1) {
    if (!requireNamespace("ABSOLUTE", quietly = TRUE)) {
        cli::cli_abort(c(
            "ABSOLUTE must be installed to use this function.",
            i = "you can install it by {.code install.packages(\"https://software.broadinstitute.org/cancer/cga/sites/default/files/data/tools/absolute/ABSOLUTE_1.0.6.tar.gz\", repos = NULL, type = \"source\")}" # nolint
        ))
    }

    # match options --------------------------------------------------------
    platform <- match.arg(
        platform, c("SNP_6.0", "Illumina_WES", "SNP_250K_STY")
    )
    copy_num_type <- match.arg(copy_num_type, c("total", "allelic"))

    if (!dir.exists(results_dir)) {
        dir.create(results_dir, recursive = TRUE)
    }

    absolute_data <- absolute_validate_seg_and_maf_data(seg = seg, maf = maf)
    absolute_filepath <- absolute_prepare_seg_and_maf_data(
        seg = absolute_data[["seg"]],
        maf = absolute_data[["maf"]],
        results_dir = results_dir
    )

    if (length(absolute_filepath[["sample_id"]]) > 0L) {
        cli::cli_inform("Running ABSOLUTE algorithm...")

        run_absolute_dir <- file.path(results_dir, "RunAbsolute")
        future.apply::future_lapply(
            absolute_filepath[["sample_id"]],
            function(sample_id) {
                maf_fn <- absolute_filepath[["maf"]][sample_id]
                if (is.null(maf_fn) || is.na(maf_fn)) {
                    maf_fn <- NULL
                }
                absolute_safe(
                    seg_dat_fn = absolute_filepath[["seg"]][sample_id],
                    maf_fn = maf_fn,
                    sample_name = sample_id,
                    sigma_p = sigma_p, max_sigma_h = max_sigma_h,
                    min_ploidy = min_ploidy, max_ploidy = max_ploidy,
                    primary_disease = primary_disease, platform = platform,
                    results_dir = run_absolute_dir,
                    max_as_seg_count = max_as_seg_count,
                    max_non_clonal = max_non_clonal,
                    max_neg_genome = max_neg_genome,
                    copy_num_type = copy_num_type,
                    min_mut_af = min_mut_af
                )
            },
            future.globals = TRUE
        )
        run_absolute_files <- file.path(
            run_absolute_dir,
            paste0(absolute_filepath[["sample_id"]], ".ABSOLUTE.RData")
        )

        run_absolute_files <- run_absolute_files[
            file.exists(run_absolute_files)
        ]

        if (!length(run_absolute_files)) {
            cli::cli_abort("No RunAbsolute results file to proceed.")
        }

        summarize_dir <- file.path(results_dir, "CreateReviewObject")

        if (dir.exists(summarize_dir)) {
            cli::cli_inform("Removing previous summarize results directory.")
            unlink(summarize_dir, recursive = TRUE)
        }
        cli::cli_inform("Summarizing multiple ABSOLUTE...")

        suppressWarnings(ABSOLUTE::CreateReviewObject(
            obj.name = "SummarizeAbsolute",
            absolute.files = run_absolute_files,
            indv.results.dir = summarize_dir,
            copy_num_type = copy_num_type,
            plot.modes = TRUE
        ))

        cli::cli_inform("Absolute auto-reviewing...")

        pp_call_fn <- file.path(
            summarize_dir,
            "SummarizeAbsolute.PP-calls_tab.txt"
        )
        modes_fn <- file.path(
            summarize_dir,
            "SummarizeAbsolute.PP-modes.data.RData"
        )

        suppressWarnings(ABSOLUTE::ExtractReviewedResults(
            reviewed.pp.calls.fn = pp_call_fn,
            analyst.id = "YJ",
            modes.fn = modes_fn,
            out.dir.base = results_dir,
            obj.name = "ReviewAbsolute",
            copy_num_type = copy_num_type,
            verbose = TRUE
        ))
        cli::cli_inform("Reviewing ABSOLUTE summary done.")
    }
}


# run_absolute utility functions --------------------------------------

absolute_safe <- function(seg_dat_fn, maf_fn,
                          sample_name, sigma_p, max_sigma_h,
                          min_ploidy, max_ploidy, primary_disease, platform,
                          results_dir, max_as_seg_count, max_non_clonal,
                          max_neg_genome, copy_num_type, min_mut_af) {
    absolute_args <- list(
        sample.name = sample_name,
        sigma.p = sigma_p, max.sigma.h = max_sigma_h,
        min.ploidy = min_ploidy, max.ploidy = max_ploidy,
        primary.disease = primary_disease, platform = platform,
        results.dir = results_dir, max.as.seg.count = max_as_seg_count,
        max.non.clonal = max_non_clonal,
        max.neg.genome = max_neg_genome,
        copy_num_type = copy_num_type,
        min.mut.af = min_mut_af
    )

    rlang::try_fetch(
        {
            suppressWarnings(rlang::inject(ABSOLUTE::RunAbsolute(
                seg.dat.fn = seg_dat_fn, maf.fn = maf_fn,
                !!!absolute_args
            )))
        },
        error = function(cnd) {
            if (any(grepl("mutations left", conditionMessage(cnd)))) {
                cli::cli_warn(c(
                    "Detecting error in sample: {.filed {sample_name}}",
                    "x" = conditionMessage(cnd),
                    "i" = "Try to fix error by removing maf ({.file maf_fn}) file"
                ))
                rlang::try_fetch(
                    {
                        suppressWarnings(rlang::inject(ABSOLUTE::RunAbsolute(
                            seg.dat.fn = seg_dat_fn, maf.fn = NULL,
                            !!!absolute_args
                        )))
                        cli::cli_inform("v" = "Fixing {.filed {sample_name}} successfully")
                    },
                    error = function(cnd2) {
                        cli::cli_warn(c(
                            "Fixing {.filed {sample_name}} failed",
                            "x" = conditionMessage(cnd2),
                            "i" = "Skipping this sample."
                        ))
                    }
                )
            } else {
                cli::cli_warn(c(
                    "Detecting error in sample: {.filed {sample_name}}",
                    "x" = conditionMessage(cnd),
                    "i" = "Skipping this sample"
                ))
            }
        }
    )
}

# validate seg and maf data to have corresponding columns -----------------
# nolint start
absolute_validate_seg_and_maf_data <- function(seg, maf = NULL) {
    if (!inherits(seg, "data.frame")) {
        cli::cli_abort(c(
            "The class of {.arg seg} must inherited from data.frame including data.frame, data.table, or tibble",
            i = "{.arg seg} is a {class(seg)}"
        ))
    }
    seg <- data.table::as.data.table(seg)

    if (!"Sample" %in% names(seg)) seg[, Sample := "SampleOne"]

    # check seg data ----------------------------------------------
    seg_cols <- c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
    if (!all(seg_cols %in% colnames(seg))) {
        cli::cli_abort(c(
            "Mising {.field columns} in {.arg maf}",
            i = "Cannot find {.field {seg_cols}}"
        ))
    }
    if (any(is.na(seg[, Sample]))) {
        cli::cli_warn(c(
            "Find NA values in {.field Sample} column of {.arg seg}",
            i = "Removing it..."
        ))
    }
    seg <- seg[!is.na(Sample), .SD, .SDcols = seg_cols]

    # check maf data ----------------------------------------------
    if (!is.null(maf)) {
        if (!inherits(maf, "data.frame")) {
            cli::cli_abort(c(
                "The class of {.arg maf} must inherited from data.frame including data.frame, data.table, or tibble",
                i = "{.arg maf} is a {class(maf)}"
            ))
        }
        maf <- data.table::as.data.table(maf)
        maf_cols <- list(
            Tumor_Sample_Barcode = "Tumor_Sample_Barcode",
            Chromosome = "Chromosome",
            Hugo_Symbol = "Hugo_Symbol",
            dbSNP_Val_Status = "dbSNP_Val_Status",
            Start_position = c("Start_position", "Start_Position"),
            t_ref_count = c("t_ref_count", "i_t_ref_count"),
            t_alt_count = c("t_alt_count", "i_t_alt_count")
        )
        idx <- vapply(maf_cols, function(x) {
            i <- match(names(maf), x, nomatch = 0L)
            if (any(i > 0L)) {
                # the first match
                if (any(i == 1L)) {
                    return(which(i == 1L, useNames = FALSE)[[1L]])
                } else {
                    return(which(i == 2L, useNames = FALSE)[[1L]])
                }
            } else {
                return(0L)
            }
        }, integer(1L))
        lack_cols <- idx == 0L
        if (any(lack_cols)) {
            cli::cli_abort(c(
                "Mising {.field columns} in {.arg maf}",
                i = "Cannot find {.field {names(maf_cols)[lack_cols]}}"
            ))
        }
        maf <- maf[, .SD, .SDcols = idx]
        data.table::setnames(maf, names(maf_cols))
    }

    lapply(list(seg = seg, maf = maf), function(x) {
        if (is.null(x)) {
            return(NULL)
        }
        x[
            , Chromosome := sub(
                pattern = "chr", replacement = "",
                as.character(Chromosome),
                perl = TRUE, ignore.case = TRUE
            )
        ]
        x[
            , Chromosome := sub(
                pattern = "X", replacement = "23",
                Chromosome, perl = TRUE, ignore.case = TRUE
            )
        ]
        x[Chromosome %in% as.character(1:23), ]
    })
}

absolute_prepare_seg_and_maf_data <- function(seg, maf = NULL, results_dir) {
    sample_id <- unique(seg[, Sample])
    if (!dir.exists(file.path(results_dir, "seg"))) {
        dir.create(file.path(results_dir, "seg"))
    }
    # prepare seg data
    seg_filepath <- file.path(results_dir, "seg", paste0(sample_id, ".seg"))
    names(seg_filepath) <- sample_id
    seg[, data.table::fwrite(
        x = .SD,
        file = seg_filepath[[unlist(.BY)]],
        sep = "\t"
    ), by = Sample]

    # prepare maf data

    if (is.null(maf)) {
        maf_filepath <- NULL
    } else {
        maf <- maf[Tumor_Sample_Barcode %in% sample_id, ]
        maf[, group_id := Tumor_Sample_Barcode]
        if (!dir.exists(file.path(results_dir, "maf"))) {
            dir.create(file.path(results_dir, "maf"))
        }
        maf_filepath <- data.table::fifelse(
            sample_id %in% maf[, Tumor_Sample_Barcode],
            file.path(results_dir, "maf", paste0(sample_id, ".maf")),
            NA_character_
        )
        names(maf_filepath) <- sample_id

        maf[, data.table::fwrite(
            x = .SD,
            file = maf_filepath[[unlist(.BY)]],
            sep = "\t"
        ), by = group_id]
    }
    list(sample_id = sample_id, seg = seg_filepath, maf = maf_filepath)
}
# nolint end
