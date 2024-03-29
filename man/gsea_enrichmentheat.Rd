% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General_functions.R
\name{gsea_enrichmentheat}
\alias{gsea_enrichmentheat}
\title{gsea_enrichmentheat}
\usage{
gsea_enrichmentheat(
  df,
  df.sourse = H1_S7.Betatraj_all.huge.clusters.GSEA,
  title,
  labelsize = 8,
  insidesize = 3,
  thedatabase = msig.db,
  headless = T,
  upsidedown = F
)
}
\arguments{
\item{df}{Selected term table}

\item{df.sourse}{The output from Fun.enrich_withFC}

\item{title}{The title of the plot}

\item{labelsize}{x/y axis label size default is 8}

\item{insidesize}{inside number label size default is 3}

\item{thedatabase}{default is msig.db}

\item{headless}{default is T}

\item{upsidedown}{default is F}
}
\value{
p
}
\description{
This function plot GSEA enrichment analysis
}
