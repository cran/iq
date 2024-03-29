\name{topN}
\alias{topN}
\title{
The topN algorithm for protein quantification
}
\description{
Estimates protein abundances using the N most intense ions.
}
\usage{
topN(X, N = 3, aggregation_in_log_space = TRUE)
}
\arguments{
\item{X}{A matrix of ion intensities in log2 space. Columns are samples and rows are fragment ions.}
\item{N}{The number of top ions used for quantification.}
\item{aggregation_in_log_space}{A logical value. If \code{FALSE}, data aggregation is performed in the original intensity space.}
}
\value{
A list of two components is returned
\item{estimate}{A vector with length equal to the number of columns of the input containing the protein abundances in log2 space.}
\item{annotation}{Reserved, currently an empty string.}
}
\references{
Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative protein
abundances from ion quantification in DIA-MS-based proteomics. \emph{Bioinformatics} 2020 Apr 15;36(8):2611-2613.
}
\author{
Thang V. Pham
}
