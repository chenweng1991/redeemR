% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VariantSummary.R
\name{add_hypermutable}
\alias{add_hypermutable}
\title{add_hypermutable
This function annotates the redeem@V.fitered with hypermutable mutations_v2 (using 3 young donor and 2 aged donors)}
\usage{
add_hypermutable(
  redeem,
  hypercut = 0.005,
  lmhc.cut1 = 40,
  lmhc.cut2 = 1.3,
  l = 583,
  r = 16190
)
}
\arguments{
\item{hypercut}{The cutoff the define the hypermutable mutations}

\item{lmhc.cut1}{Also annotate with the "LMHC", lmhc.cut1 cutoff by number of cells sharing a mutation}

\item{lmhc.cut2}{Also annotate with the "LMHC", lmhc.cut2 cutoff by mean count}

\item{redeemR}{}
}
\description{
add_hypermutable
This function annotates the redeem@V.fitered with hypermutable mutations_v2 (using 3 young donor and 2 aged donors)
}
