% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuidTree.R
\name{FromDist2Graph}
\alias{FromDist2Graph}
\title{FromDist2Graph
From distance object or matrix to graph, default is to return igraph object
This function was developed based on}
\usage{
FromDist2Graph(d, k.param = 30, return_igraph = T)
}
\arguments{
\item{d}{the distance matrix,  this can be either dist or a matrix}

\item{k.param}{K default is 30}

\item{return_igraph}{Wheather return igraph, default is T which return igraph. Otherwise, return adjacent matrix}
}
\value{
igraph or adjacent matrix
}
\description{
FromDist2Graph
From distance object or matrix to graph, default is to return igraph object
This function was developed based on
}
