% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_.R
\name{sample_sim}
\alias{sample_sim}
\title{Sample features with similar expression level}
\usage{
sample_sim(matrix, features, size = 10^3L, nbin = 25L)
}
\arguments{
\item{matrix}{Normalized gene expression matrix with row is genes and column
is samples.}

\item{features}{The template for gene sampling, considering genes with
comparable expression levels of these genes.}

\item{size}{A non-negative integer giving the number of sampling.}

\item{nbin}{Divides the dataset into approximately \code{nbin} clusters based on
gene expression levels. Subsequently, the sampling process will be applied
within each cluster for the specified features.}
}
\value{
A list of character for each sampled genes.
}
\description{
Sample features with similar expression level
}
\references{
Barkley, D., Moncada, R., Pour, M. et al. Cancer cell states recur across
tumor types and form specific interactions with the tumor microenvironment.
Nat Genet 54, 1192–1201 (2022). https://doi.org/10.1038/s41588-022-01141-9
}
\seealso{
\url{https://github.com/yanailab/PanCancer}
}
