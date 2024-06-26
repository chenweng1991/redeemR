% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuidTree.R
\name{quick_w_cosine}
\alias{quick_w_cosine}
\title{Compute weighted cosine distance}
\usage{
quick_w_cosine(M, w)
}
\arguments{
\item{M}{the binary matrix, Each row is a cell, each column is a variant, generated by Make_matrix}

\item{w}{weight for each variant, a vector}
}
\value{
dist object
}
\description{
Compute weighted cosine distance
}
