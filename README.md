
# scMitoTracing

<!-- badges: start -->
<!-- badges: end -->

The goal of scMitoTracing is to analyze the deep mito variants and lineage tracing
Package testing Jupyter notebook is here lab/solexa_weissman/cweng/Packages/scMitoTracing/scMitoTracingRunningTest.ipynb

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
DN1CD34_1.depth<-DepthSummary(WD)
DN1CD34_1.VariantsGTSummary<-CW_mgatk.read(WD,Processed =T)
DN1CD34_1.Variants.feature.lst<-Vfilter_v3(InputSummary=DN1CD34_1.VariantsGTSummary,depth=DN1CD34_1.depth)
```

### Plotting for variants

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

## Advanced options for variant parsing

### DepthSummary from a subset of cells
```r
library(scMitoTracing)
library(EZsinglecell2)
WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/MTenrichCombine/Enrich/CW_mgatk/final"
## Read in the CD34 Seurat object
Donor01_CD34_1_Multiome_wrapper<-readRDS("/lab/solexa_weissman/cweng/Projects/RDS/Donor01_CD34_1_Multiome_wrapper.RDS")
## Translate the RNA names into ATAC/Mito names
Donor01_CD34_1_Multiome.meta<-Translate_RNA2ATAC(Donor01_CD34_1_Multiome_wrapper$seurat@meta.data)
## Select HSC cells and parse the HSC mito depth information
HSCcells<-subset(Donor01_CD34_1_Multiome.meta,CellType=="HSC")$ATACName
DN1CD34_1.depth<-DepthSummary(WD,Processed = F,CellSubset = HSCcells,cellSubSetName = "HSC")
```

## Making tree
### Below is an example that combine both CD34 and BMMC data into one mitoTracing object

``` r
library(EZsinglecell2)
## CD34 This is the same as above
WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/MTenrichCombine/Enrich/CW_mgatk/final"
#DN1CD34_1.depth<-DepthSummary(WD)
DN1CD34_1.depth<-readRDS(paste(WD,"/DN1CD34_1.depth.RDS",sep=""))
DN1CD34_1.VariantsGTSummary<-CW_mgatk.read(WD,Processed =T)
DN1CD34_1.Variants.feature.lst<-Vfilter_v3(InputSummary=DN1CD34_1.VariantsGTSummary,depth=DN1CD34_1.depth)

## BMMC1
WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/Run210706L2/COMBINE_CALLV/CW_mgatk_Full/final"
# DN1BMMC1.depth<-DepthSummary(WD)
#saveRDS(DN1BMMC1.depth,paste(WD,"/DN1BMMC1.depth.RDS",sep=""))
DN1BMMC1.depth<-readRDS(paste(WD,"/DN1BMMC1.depth.RDS",sep=""))
DN1BMMC1.VariantsGTSummary<-CW_mgatk.read(WD,Processed =T)
DN1BMMC1.Variants.feature.lst<-Vfilter_v3(InputSummary=DN1BMMC1.VariantsGTSummary,depth=DN1BMMC1.depth)

## Make meta
bmmc.filtered<-readRDS("/lab/solexa_weissman/cweng/Projects/RDS/mutiome.donor01_bmmc01.RNA.filtered.RDS") 
Donor01_CD34_1_Multiome_wrapper<-readRDS("/lab/solexa_weissman/cweng/Projects/RDS/Donor01_CD34_1_Multiome_wrapper.RDS")
DN1BMMC1.meta<-Translate_RNA2ATAC(meta=bmmc.filtered@meta.data[,c("nCount_RNA","nFeature_RNA","CellType")])
DN1CD34_1.meta<-Translate_RNA2ATAC(meta=Donor01_CD34_1_Multiome_wrapper$seurat@meta.data[,c("nCount_RNA","nFeature_RNA","CellType")])

## Make CD34_BMMC_mitoTracing.VerySensitive
CD34_BMMC_mitoTracing.VerySensitive<-Create_mitoTracing(GTsummary_list=list(DN1CD34_1.VariantsGTSummary,DN1BMMC1.VariantsGTSummary),
                depth_list=list(DN1CD34_1.depth,DN1BMMC1.depth),
                feature.list_list=list(DN1CD34_1.Variants.feature.lst,DN1BMMC1.Variants.feature.lst),
                meta_list=list(DN1CD34_1.meta,DN1BMMC1.meta),
                thr="VerySensitive",
                labels=c("CD34","BMMC"))

## Show basics of the CD34_BMMC_mitoTracing.VerySensitive
CD34_BMMC_mitoTracing.VerySensitive
```

### Make variant matrix and clonal clustering
``` r
CD34_BMMC_mitoTracing.VerySensitive<-Make_matrix(CD34_BMMC_mitoTracing.VerySensitive)
CD34_BMMC_mitoTracing.VerySensitive<-SeuratLSIClustering(CD34_BMMC_mitoTracing.VerySensitive)
CD34_BMMC_mitoTracing.VerySensitive<-AddDatatoplot_clustering(CD34_BMMC_mitoTracing.VerySensitive)
```


