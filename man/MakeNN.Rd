% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuidTree.R
\name{MakeNN}
\alias{MakeNN}
\title{Define a function to make nn list, which can be further used to make adjacency matrix
This scan row by row, looking for k.param nearest neighbours}
\usage{
MakeNN(d, k.param = 15)
}
\arguments{
\item{d}{Distance matrix, can be a dist object or matrix}

\item{k.param}{Default is 15}
}
\value{
return an nn list, which has two components: nn$idx and nn$dist
}
\description{
Define a function to make nn list, which can be further used to make adjacency matrix
This scan row by row, looking for k.param nearest neighbours
}
