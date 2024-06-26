% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuidTree.R
\name{Run_Lin_regression_poi}
\alias{Run_Lin_regression_poi}
\title{Run_Lin_regression_poi
Firstly used in HSC_multiome_Het_2.ipynb
This function was developed based on}
\usage{
Run_Lin_regression_poi(
  LinOut,
  regress_factor = c("OutLevel.scale", "OutLevel_NPadj.scale", "Lym", "Mye", "MK", "ME"),
  n.cores = 8
)
}
\arguments{
\item{LinOut}{produced by MakeDF4Regress}

\item{regress_factor}{default is c("OutLevel.scale","OutLevel_NPadj.scale","Lym","Mye","MK","ME")}

\item{n.cores}{=8}
}
\description{
Run_Lin_regression_poi
Firstly used in HSC_multiome_Het_2.ipynb
This function was developed based on
}
