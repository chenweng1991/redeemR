% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuidTree.R
\name{Subset_redeemR}
\alias{Subset_redeemR}
\title{Subset_redeemR Subset a redeemR object by selecting a subset of cells, return a new redeemR object with only 4 slots: para; CellMeta; Cts.Mtx.bi; UniqueV, can be used for downstreme compute distance, clonal clustering, make tree, etc}
\usage{
Subset_redeemR(redeemR, Cells, ExtraInfo = "Subset from ... ")
}
\arguments{
\item{redeemR}{The Parent redeemR object eg. DN4_HSC_redeemR.Sensitive}

\item{Cells}{Important, give a vector of Cell names(ATAC cell names)}

\item{ExtraInfo}{Extra information, usually "Subset from ..."}
}
\value{
redeemR Object
}
\description{
Subset_redeemR Subset a redeemR object by selecting a subset of cells, return a new redeemR object with only 4 slots: para; CellMeta; Cts.Mtx.bi; UniqueV, can be used for downstreme compute distance, clonal clustering, make tree, etc
}
