
# scMitoTracing

<!-- badges: start -->
<!-- badges: end -->

The goal of scMitoTracing is to ...

## Installation

You can install the development version of scMitoTracing like so:

``` r
devtools::install_github("chenweng1991/scMitoTracing")
library(EZsinglecell2)
```

## Quick start


### Variant parsing and analysis

``` r
library(scMitoTracing)
WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/MTenrichCombine/Enrich/CW_mgatk/final"
DN1CD34_1.depth<-DepthSummary(WD,Processed = T)
DN1CD34_1.VariantsGTSummary<-CW_mgatk.read(WD,Processed =T)
DN1CD34_1.Variants.feature.lst<-Vfilter_v3(InputSummary=DN1CD34_1.VariantsGTSummary,depth=DN1CD34_1.depth)
```

### Plotting

Plot for depth
``` r
plot_depth(DN1CD34_1.depth$Total,"Total")
plot_depth(DN1CD34_1.depth$VerySensitive,"VerySensitive")
```

Plot for variant signature
``` r
options(repr.plot.width=16, repr.plot.height=6)
ps<-list()
for(name in names(DN1CD34_1.Variants.feature.lst)){
p<-MutationProfile.bulk(DN1CD34_1.Variants.feature.lst[[name]]$Variants)+ggtitle(name)+theme(title =element_text(size=20))
ps<-c(ps,list(p))
}
grid.arrange(grobs=ps)
```

Plot rssential variant metrics
```r
plot_variant(DN1CD34_1.VariantsGTSummary,DN1CD34_1.Variants.feature.lst,depth=DN1CD34_1.depth,cat=c("Total","VerySensitive","Sensitive","Specific"),p4xlim = 30)
```
