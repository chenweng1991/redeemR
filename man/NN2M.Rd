% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuidTree.R
\name{NN2M}
\alias{NN2M}
\title{Define a function convert nn list to adjacency matrix that can be further used for igraph}
\usage{
NN2M(nn)
}
\arguments{
\item{nn}{nn list, which has two components: nn$idx and nn$dist}
}
\value{
return an nn.matrix. This is adjacency matrix can be input to igraph graph<-graph_from_adjacency_matrix(nn.matrix,diag = F,mode = "undirected")
}
\description{
Define a function convert nn list to adjacency matrix that can be further used for igraph
}
