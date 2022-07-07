% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VariantSummary.R
\name{ComputeRejectRate}
\alias{ComputeRejectRate}
\title{Function to compute the reject rate(The filtering rate in concensus variant calling)}
\usage{
ComputeRejectRate(WD)
}
\arguments{
\item{WD}{The path to the work space usually  XXX/mitoV/final}
}
\value{
a dataframe that store the percentage of variant in a given threahold again total
}
\description{
This function allows you to computae the filtering rate for each single cell
}
\examples{
DN9_BMMC_RejectRate<-ComputeRejectRate("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9/Donor9/DN9_BMMC/MTenrichCombine/mitoV/final/")
}