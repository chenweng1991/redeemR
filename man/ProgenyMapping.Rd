% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BuidTree.R
\name{ProgenyMapping}
\alias{ProgenyMapping}
\title{Define a function to perform single-cell based hard porogeny assignment
This function was developed based on DN4T2.basics.ipynb}
\usage{
ProgenyMapping(
  HSC_redeemR = DN4_PhenoHSC_redeemR.verysensitive,
  Full_redeemR = DN4_BMMC_HSPC_HSC_redeemR.verysensitive,
  distCut = 0.95,
  d = "w_jaccard"
)
}
\arguments{
\item{HSC_redeemR}{The HSC_redeemR is the redeemR object for defined HSC}

\item{Full_redeemR}{The FULL_redeemR is the redeemR object for the full BMMC_HSPC_HSC}

\item{distCut}{Default is 0.95, the distance, below which I define as the related progeny}
}
\description{
Define a function to perform single-cell based hard porogeny assignment
This function was developed based on DN4T2.basics.ipynb
}
