![GitHub R package version](https://img.shields.io/github/r-package/v/chenweng1991/redeemR?label=ReDeeM)
![GitHub last commit](https://img.shields.io/github/last-commit/chenweng1991/redeemR)
![GitHub](https://img.shields.io/github/license/chenweng1991/redeemR)

# redeemR

<!-- badges: start -->
<!-- badges: end -->


ReDeeM: single-cell **Re**gulatory multi-omics with **Dee**p **M**itochondrial mutation profiling. ReDeeM is an engineering-free single cell sequencing technology for somatic mutation-based retrospective lineage tracing, clonal/subclonal mapping, with simultaneous cell-state profiling (transcriptomics and chromatin accessibility).


The redeemR is the second component in ReDeeM computational toolsets
- [redeemV](https://github.com/chenweng1991/REDEEM-V) is an in-house python package for mapping and deep mitochondrial mutation calling. (Input from rawdata)
- [redeemR](https://github.com/chenweng1991/REDEEM-R) is an in-house R package for downstream lineage tracing and single cell integrative analysis. (Input from redeemV)

The goal of redeemR is to refine the somatic mitochondiral mutations, build the single-cell phylogenetic tree and facilitate clonal/subclonal-resolved single-cell multiomics analysis.


## Installation

You can install the development version of redeemR:

``` r
devtools::install_github("chenweng1991/redeemR")
library(redeemR)
```

## Documentation 
- [Getting started]


## Contact
If you run into issues and would like to report them, please submit an Issue. Also please feel free to contact authors: cweng{at}broadinstitute.org, fyu{at}broadinstitute.org.
