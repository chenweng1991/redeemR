% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Plot.R
\name{plot_variant_legacy}
\alias{plot_variant_legacy}
\title{Legacy Function to plot variant metrics}
\usage{
plot_variant_legacy(
  GTSummary,
  feature.list,
  depth,
  cat = c("Total", "VerySensitive", "Sensitive", "Specific"),
  p4xlim = 50,
  QualifyCellCut = 10
)
}
\arguments{
\item{GTSummary}{GTSummary from CW_mgatk.read}

\item{feature.list}{feature.list from Vfilter_v3}

\item{cat}{The catogories, it can be cat = c("Total", "VerySensitive", "Sensitive", "Specific") or a subset}

\item{p4xlim}{the p4 xlim(number of variant per cell), default is 50}

\item{QualifyCellCut}{median coverage for qualified cells, default is 10}
}
\value{
no returns, directly plot
}
\description{
This function works with CW_mgatk.read and Vfilter_v3
This allows you to plot the mito mutation metrics
This legacy function is useful to look at all threadhold simultaneous
For each category(stringency),
p1: Variant allele frequency(VAF);
p2: Heteroplasmy histogram
p3: CellN(Number of caells that carry the variants) VS maxcts( The number of variant counts in the highest cell)
p4: Histogram to show the distribution of the number of variant per cell
}
\examples{
plot_variant(DN1CD34_1.VariantsGTSummary,DN1CD34_1.Variants.feature.lst,depth=DN1CD34_1.depth,cat=c("Total","VerySensitive","Sensitive","Specific"),p4xlim = 30)
}
