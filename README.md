
# redeemR

<!-- badges: start -->
<!-- badges: end -->

The use of redeemR is to analyze the deep mito variants and lineage tracing
Package testing Jupyter notebook is here 

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
### Clonal clustering and distance metrics (Updated 3/30/22)
``` r
## Make a weight matrix (For different variants) based on recurrencey. It is required for weighted distance calculation
data(CellPCT)
V.weight<-data.frame(weight=1-CellPCT$muRate)
V.weight$Variants<-paste("Variants",gsub("_","",CellPCT$Variant),sep="")

BMMC_mitoTracing<-SeuratLSIClustering(BMMC_mitoTracing,lsidim=2:50,rmvariants=c("Variants310TC","Variants3109TC","Variants5764CT"))
BMMC_mitoTracing<-AddDist(BMMC_mitoTracing,weightDF=V.weight,LSIdist=T,dim=2:50)
```

### Subset a mitotracing object (Updated 5/8/22)
``` r
# Example from HSC_multiome_Het.ipynb
DN4_stemcell_mitoTracing.seed.verysensitive<-Subset_MitoTracing(DN4_stemcell_mitoTracing.verysensitive,Cells=topcell)
DN4_stemcell_mitoTracing.seed.verysensitive<SeuratLSIClustering(DN4_stemcell_mitoTracing.seed.verysensitive,res=2,lsidim=3:50,rmvariants=c("Variants310TC","Variants3109TC","Variants5764CT"))
DN4_stemcell_mitoTracing.seed.verysensitive<-AddDist(DN4_stemcell_mitoTracing.seed.verysensitive,weightDF=V.weight,LSIdist=T,dim=3:50)
```

### Make tree (Updated 5/8/22)
```r
DN4_stemcell_mitoTracing.seed.verysensitive<-Make_tree(DN4_stemcell_mitoTracing.seed.verysensitive,d = "w_jaccard",algorithm = "nj",onlyreturntree = F)
```

### Add depth, assign variant to tree (Updated 5/8/22)  
``` r
## Load depth
## Prepare the depth, given the order created: "BMMC":-1 ,"HSPC": -2, "HSC": -3
#DN4_WD_BMMC="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9/Donor4/DN4_BMMC/MTenrichCombine/mitoV/final/"
DN4_WD_HSPC="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9/Donor4/DN4_HPC/MTenrichCombine/mitoV/final/"
DN4_WD_HSC="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9/Donor4/DN4_HSC/MTenrichCombine/mitoV/final/"
q_HSPC<-read.table(paste(DN4_WD_HSPC,"QualifiedTotalCts",sep=""))
q_HSC<-read.table(paste(DN4_WD_HSC,"QualifiedTotalCts",sep=""))
q_HSPC$V1<-paste(q_HSPC$V1,"2",sep="_")
q_HSC$V1<-paste(q_HSC$V1,"3",sep="_")
quality<-rbind(q_HSPC,q_HSC)
DN4_stemcell.seed.depth<-subset(quality,V1 %in% DN4_stemcell_mitoTracing.seed.verysensitive@CellMeta$Cell)

## add depth
DN4_stemcell_mitoTracing.seed.verysensitive<-Add_DepthMatrix(DN4_stemcell_mitoTracing.seed.verysensitive,quality)

## maximum-liklihood-based variant assignment
DN4_stemcell_mitoTracing.seed.verysensitive<-Add_AssignVariant(DN4_stemcell_mitoTracing.seed.verysensitive,n.cores = 8)
```

### Make Allnodes(Node|Parent|Freq|CladeSize), visulize the branchlength by number of assigned variants (Updated 5/8/22)
``` r
Allnodes=MakeAllNodes(mitoTracing,prob.cut=0.1)
TreeData<-DN4_stemcell_mitoTracing.seed.verysensitive@TREE@phylo
TreeData$edge.length<-Allnodes[match(DN4_stemcell_mitoTracing.seed.verysensitive@TREE@phylo$edge[,2],Allnodes$Node),] %>% .$Freq
options(repr.plot.width=3, repr.plot.height=24,repr.plot.res=120)
as.treedata(TreeData) %>% ggtree(.)+theme_tree2()
```

### Cut the tree in to clones(Updated 5/8/22)  
A useful tool from phanforn[https://rdrr.io/cran/phangorn/man/Ancestors.html] 

``` r
DN4_stemcell_mitoTracing.seed.verysensitive<-Add_tree_cut(DN4_stemcell_mitoTracing.seed.verysensitive,MinCell = 30,N = 1, prob.cut = 0.3)

#Visulize
DN4_stemcell_mitoTracing.seed.verysensitive@CellMeta$Clone_merge<-as.factor(DN4_stemcell_mitoTracing.seed.verysensitive@CellMeta$Clone_merge)
library(RColorBrewer)
n=length(unique(DN4_stemcell_mitoTracing.seed.verysensitive@CellMeta$Clone_merge))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
options(repr.plot.width=6, repr.plot.height=6,repr.plot.res=100)
p<-ggtree(DN4_stemcell_mitoTracing.seed.verysensitive@TREE@treedata,layout="circular", branch.length='none') 
p+geom_fruit( 
         data=DN4_stemcell_mitoTracing.seed.verysensitive@CellMeta, 
         geom=geom_tile, 
         mapping=aes(y=Cell,x=2,fill=Clone_merge), 
         pwidth=0.001, 
         width=3, 
         offset=0.05
     )+scale_fill_manual(values=col_vector)
```



