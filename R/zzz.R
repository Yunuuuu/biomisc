.onAttach <- function(libname, pkgname) {
    version <- utils::packageDescription(pkgname, fields = "Version")

    packageStartupMessage(
        "Package: ", pkgname, " (version: ", version, ")\n",
        "Miscellaneous functions for Bioinformatics!"
    )
    invisible()
}

if (getRversion() >= "2.15.1") {
    utils::globalVariables(
        c( # run_absolute
            "Chromosome", "Sample", "Tumor_Sample_Barcode",
            "group_id"
        )
    )
}
