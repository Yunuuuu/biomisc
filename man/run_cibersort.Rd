% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_cibersort.R
\name{run_cibersort}
\alias{run_cibersort}
\title{Cell-type identification by estimating relative subsets of RNA transcripts
(CIBERSORT)}
\usage{
run_cibersort(
  mixture_data,
  sig_data = NULL,
  perm = 200L,
  quantile_norm = TRUE,
  absolute = FALSE,
  abs_method = "sig_score"
)
}
\arguments{
\item{mixture_data}{a \code{matrix} object whose rownames is the gene names and is
consistent with \code{sig_data} gene names; Data should be in non-log space.
Note: if maximum expression value is <50; CIBERSORT will assume that data are
in log space, and will anti-log all expression values by 2x. If gene symbols
are redundant, CIBERSORT will choose the one with highest mean expression
across the mixtures. CIBERSORT performs a feature selection and therefore
typically does not use all genes in the signature matrix. It is generally ok
if some genes are missing from the user’s mixture file. If <50\% of signature
matrix genes overlap, CIBERSORT will issue a warning.  Normal quantification
of RNA-seq like FPKM, and TPM can be used.}

\item{sig_data}{CIBERSORT requires an input matrix of reference gene
expression signatures, or signature matrix with gene names in rownames, for
routine analysis. This is stored in a Signature Genes File and consists of a
table with groups of "barcode" genes whose expression values collectively
define a unique gene expression signature for each component pure cell
population that will be used to deconvolute the mixture. if \code{NULL},
\code{LM22} with HUGO gene symbols as gene names will be used.}

\item{perm}{Number of permutations; set to >=100 to calculate p-values.
default: \code{200}.}

\item{quantile_norm}{Quantile normalization of input mixture. default: \code{TRUE}}

\item{absolute}{Run CIBERSORT in absolute mode default: \code{FALSE}. note
that cell subsets will be scaled by their absolute levels and will not be
represented as fractions (to derive the default output, normalize absolute
levels such that they sum to 1 for each mixture sample); the sum of all
cell subsets in each mixture sample will be added to the ouput ('Absolute
score'). If LM22 is used, this score will capture total immune content.}

\item{abs_method}{if absolute is set to \code{TRUE}, abs_method choose
method: 'no_sumto1' or 'sig_score' - sig_score = for each mixture sample,
define S as the median expression level of all genes in the signature
matrix divided by the median expression level of all genes in the mixture.
Multiple cell subset fractions by S. - no_sumto1 = remove sum to 1
constraint}
}
\value{
A \link[data.table:data.table]{data.table}
}
\description{
Robust enumeration of cell subsets from tissue expression
profiles
}
\details{
By default, CIBERSORT estimates the relative fraction of each cell
type in the signature matrix, such that the sum of all fractions is equal
to 1 for a given mixture sample.

CIBERSORT can also be used to produce a score that quantitatively measures
the overall abundance of each cell type (as described in
\href{https://www.nature.com/articles/nmeth.3337}{Analysis of deconvolution
consistency}). Briefly, the absolute immune fraction score is estimated by
the median expression level of all genes in the signature matrix divided by
the median expression level of all genes in the mixture. more details see
\href{https://cibersort.stanford.edu/manual.php#run}{CIBERSORT Manual}
}
\references{
\itemize{
\item Newman, A., Liu, C., Green, M. et al. Robust enumeration of cell subsets
from tissue expression profiles. Nat Methods 12, 453–457 (2015).
\url{https://doi.org/10.1038/nmeth.3337}
\item Chen B., Khodadoust M.S., Liu C.L., Newman A.M., Alizadeh A.A. (2018)
Profiling Tumor Infiltrating Immune Cells with CIBERSORT. In: von Stechow
L. (eds) Cancer Systems Biology. Methods in Molecular Biology, vol 1711.
Humana Press, New York, NY.  \url{https://doi.org/10.1007/978-1-4939-7493-1_12}
}
}
