% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuidTree.R
\name{ProgenyMapping_np}
\alias{ProgenyMapping_np}
\title{ProgenyMapping_np
Define a function to compute network propagation based probability
FromDist2Graph is needed to convert fistance matrix into MNN graph}
\usage{
ProgenyMapping_np(
  HSC_mitoTracing = DN4_stemcell_mitoTracing.seed.verysensitive,
  Full_mitoTracing = DN4_BMMC_HSPC_HSC_mitoTracing.verysensitive,
  CloneCol = "Clone_merge",
  k = 30,
  gm = 0.5,
  useLSI = F,
  useSCAVENGE_LSI = F,
  subsample = F,
  ProbCut = 0.7,
  Celltype = "Rig.CellType"
)
}
\arguments{
\item{HSC_mitoTracing}{The HSC_mitoTracing is the mitoTracing object for defined HSC,  have already gone through Add_DepthMatrix--Add_AssignVariant--Add_tree_cut, otherwise, need
othereise, need a column in CellMeta that indicates the clone ID}

\item{Full_mitoTracing}{The FULL_mitoTracing is the mitoTracing object for the full BMMC_HSPC_HSC}

\item{CloneCol}{"Clone_merge"}

\item{k}{the k.param used for MNN graph}

\item{gm}{gamma default is 0.05 which mean 95\% information is passing out}

\item{ProbCut}{The cutoff of the maximum probability for a given progeny cell(If the maximum probability is lower than ProbCut, it will be filtered)}

\item{Celltype}{The column to be used in aggregate into lineages}
}
\value{
a list of two  ALLmeta.npClone (A meta data with last column npClone), np_mat (the network propagation matrix))
}
\description{
ProgenyMapping_np
Define a function to compute network propagation based probability
FromDist2Graph is needed to convert fistance matrix into MNN graph
}