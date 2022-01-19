#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' An intermediate  S4 class Datatoplots
#'
#' @slot clustering dataframe that store the data to plot
#' @export
Datatoplots<-setClass(
    "Datatoplots",
    slots=c(
    clustering="data.frame"
    )
)


#' An intermediate  S4 class Datatoplots
#'
#' @slot jaccard distance object dist: Jaccard distance
#' @slot Dice distance object dist: Dice distance
#' @slot jaccard3W distance object dist: jaccard3W
#' @import stats
#' @export
DistObjects<-setClass(
    "DistObjects",
    slots=c(
    jaccard="dist",
    Dice="dist",
    jaccard3W="dist"
    )
)


#' An intermediate  S4 class Tree that store tree info
#'
#' @slot phylo the phylo tree class from ape package
#' @slot treedata treedata class from tidytree
#' @slot records character to store annotations
#' @import ape phytools phangorn treeio ggtree tidytree ggtreeExtra
#' @export
TREE<-setClass(
    "TREE",
    slots=c(
    phylo="phylo",
    treedata="treedata",
    records="character"
    )
)


#' Major mitoTracing class that store clonal-resolved multi-omics
#'
#' @slot GTsummary.filtered  The Mitochondrial genotype data frame
#' @slot CellMeta Store meta data for each cell type
#' @slot V.fitered.list a list of data frame of variant metrics, VAF, cellN, etc (each for different stringency),
#' @slot UniqueV A character showing the number of usable variant
#' @slot Cts.Mtx A sparse matrix cell-mitoVariants, store the variant count
#' @slot Cts.Mtx.bi A sparse matrix cell-mitoVariants, The variant count has been binarized into 0 and 1
#' @slot para A character showing the parameter of this object
#' @slot Seurat Seurat object storing the clonal clustering results
#' @slot DataToplotList The customized class of Datatoplots: A list of dataframe for further plotting
#' @slot DistObjects The customized class that stores the cell-cell distances
#' @slot TREE The customized class that wraps phylogenetic tree
#' @import Seurat
#' @export
mitoTracing<-setClass(
    "mitoTracing",
    slots=c(GTsummary.filtered="data.frame",
            CellMeta="data.frame",
            V.fitered.list="list",
            UniqueV="character",
            Cts.Mtx="dgCMatrix",
            Cts.Mtx.bi="dgCMatrix",
            para="character",
            Seurat="Seurat",
            DataToplotList="Datatoplots",
            DistObjects="DistObjects",
            TREE="TREE"
           )
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Make_matrix
#' This will make the matixies of Cell VS mitochondrial variants and return mitoTracing
#' Results stored in Cts.Mtx and Cts.Mtx.bi
#' @param object mitoTracin class
#' @export
setGeneric(name="Make_matrix", def=function(object)standardGeneric("Make_matrix"))

#' SeuratLSIClustering
#' This will use the mito variants for Seurat clustering (LSI based)
#' @param object mitoTracin class
#' @export
setGeneric(name="SeuratLSIClustering", def=function(object,...)standardGeneric("SeuratLSIClustering"))

#' AddDatatoplot_clustering
#' This prepare the clonal clustering data to plot
#' @param object mitoTracin class
#' @export
setGeneric(name="AddDatatoplot_clustering", def=function(object,...)standardGeneric("AddDatatoplot_clustering"))

#' AddDist
#' This add Jaccard, Dice, Jaccard3W distance and stored in DistObjects
#' @param object mitoTracin class
#' @export
setGeneric(name="AddDist", def=function(object,...) standardGeneric("AddDist"))

#' Make_tree
#' This will generate a basic phylogenetic tree
#' @param object mitoTracin class
#' @param d "jaccard" or "Dice" or "jaccard3W"
#' @param algorithm the algorithm used to build the tree, choose from "nj" and "upgma"
#' @export
setGeneric(name="Make_tree", def=function(object,d="jaccard", algorithm="upgma",onlyreturntree=F,...) standardGeneric("Make_tree"))

#' Add_Tree
#' Optional, if a phylogentic tree object phylo is already available, can be directly added to the mitoTracing
#' @param object mitoTracin class
#' @param phylo phyogenetic tree object
#' @export
setGeneric(name="AddTree", def=function(object,phylo,...) standardGeneric("AddTree"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' show
#' This will show the basics of mitoTracin class
#' @param object mitoTracin class
#' @return print out basics
#' @export
setMethod(f="show",
          signature="mitoTracing",
          definition=function(object){
              print(object@para)
              print(paste("Total Cell number:",nrow(object@CellMeta)))
              print(table(object@CellMeta$Label))
              print(paste("Total Variant number:",length(object@UniqueV)))
              print(paste("Slot:",slotNames(object)))
          })

#' Make_matrix
#' This will make the matixies of Cell VS mitochondrial variants and return mitoTracing
#' Results stored in Cts.Mtx and Cts.Mtx.bi
#' @param object mitoTracin class
#' @return mitoTracin class
#' @export
setMethod(f="Make_matrix",
          signature="mitoTracing",
          definition=function(object){
                require(dplyr)
                require(Matrix.utils)
                Cts.Mtx<-dMcast(object@GTsummary.filtered,Cell~Variants,value.var = "Freq")
                colnames(Cts.Mtx)<-strsplit(as.character(colnames(Cts.Mtx)),"_") %>% sapply(.,function(x){paste(x[1],x[2],x[3],sep="")})
                Cts.Mtx.bi<-Cts.Mtx
                Cts.Mtx.bi[Cts.Mtx.bi>=1]<-1
                object@Cts.Mtx.bi<-Cts.Mtx.bi
                object@Cts.Mtx<-Cts.Mtx
#                 validObject(object)
                return(object)
})

#' Make_tree
#' This will generate a basic phylogenetic tree
#' @param object mitoTracin class
#' @param d "jaccard" or "Dice" or "jaccard3W"
#' @param algorithm the algorithm used to build the tree, choose from "nj" and "upgma"
#' @return mitoTracin class
#' @export
setMethod(f="Make_tree",
          signature="mitoTracing",
          definition=function(object,d,algorithm,onlyreturntree=F){
          dist<-slot(object@DistObjects,d)
          if(algorithm=="nj"){
          phylo<-nj(dist)
          }else if (algorithm=="upgma"){
          phylo<-upgma(dist)
          }
          treedata<-as.treedata(phylo)
          TREEobject<-new("TREE",phylo=phylo,treedata=treedata,records=paste(d,algorithm,sep="-"))
          if(onlyreturntree){
          return(TREEobject)
          }else{
          object@TREE<-TREEobject
          return(object)
          }
})


#' SeuratLSIClustering
#' This will use the mito variants for Seurat clustering (LSI based)
#' @param  mitoTracing class
#' @param binary  Default is tree, to make use of the binary matrix
#' @param res     Default os 0.3, the resolution of the clustering
#' @return mitoTracing class
#' @export
setMethod(f="SeuratLSIClustering",
          signature="mitoTracing",
          definition=function(object,binary=T,res=0.3){
          require(Signac)
          if(binary){
              Cell_Variant.seurat<-CreateSeuratObject(counts = t(as.matrix(object@Cts.Mtx.bi)), assay = "mitoV")
          }else{
              Cell_Variant.seurat<-CreateSeuratObject(counts = t(as.matrix(object@Cts.Mtx)), assay = "mitoV")
          }
          VariableFeatures(Cell_Variant.seurat) <- row.names(Cell_Variant.seurat) #names(which(Matrix::rowSums(Cell_Variant.seurat) > 100))
          Cell_Variant.seurat <- RunTFIDF(Cell_Variant.seurat, n = 50)
          Cell_Variant.seurat<- FindTopFeatures(Cell_Variant.seurat, min.cutoff = 'q0')
          Cell_Variant.seurat <- RunSVD(Cell_Variant.seurat, n = 50)
          Cell_Variant.seurat <- RunUMAP(Cell_Variant.seurat, reduction = "lsi", dims = 1:20)
          Cell_Variant.seurat <- FindNeighbors(Cell_Variant.seurat,reduction ="lsi"  ,dims = 1:20)
          Cell_Variant.seurat <- FindClusters(Cell_Variant.seurat, resolution = 0.3)
          object@Seurat<-Cell_Variant.seurat
          return(object)
})

#' AddDatatoplot_clustering
#' This prepare the clonal clustering data to plot
#' @param object mitoTracin class
#' @return mitoTracing class
#' @export
setMethod(f="AddDatatoplot_clustering",
          signature="mitoTracing",
          definition=function(object){
          row.names(object@CellMeta)<-object@CellMeta$Cell
          datatoplot<-Tomerge_v2(object@Seurat@meta.data,object@Seurat@reductions$umap@cell.embeddings) %>% Tomerge_v2(.,object@Seurat@reductions$lsi@cell.embeddings[,1:6]) %>% Tomerge_v2(.,object@CellMeta)
          object@DataToplotList<-Datatoplots()
          object@DataToplotList@clustering<-datatoplot
          return(object)
})

#' AddDist
#' This add Jaccard, Dice, Jaccard3W distance and stored in DistObjects
#' @param object mitoTracin class
#' @return mitoTracing class
#' @export
setMethod(f="AddDist",
          signature="mitoTracing",
          definition=function(object){
          d.Jaccard<-BinaryDist(object@Cts.Mtx.bi,method="Jaccard")
          d.Dice<-BinaryDist(object@Cts.Mtx.bi,method="Dice")
          d.3WJaccard<-BinaryDist(object@Cts.Mtx.bi,method="3WJaccard")
          object@DistObjects<-new("DistObjects",jaccard=d.Jaccard, Dice=d.Dice,jaccard3W=d.3WJaccard)
          return(object)
          })

#' Add_Tree
#' Optional, if a phylogentic tree object phylo is already available, can be directly added to the mitoTracing class in slot TREE
#' @param object mitoTracin class
#' @param phylo phyogenetic tree object
#' @param object mitoTracin class
#' @return mitoTracing class
#' @export
setMethod(f="AddTree",
          signature="mitoTracing",
          definition=function(object,phylo,record=""){
          TREEobject<-new("TREE",phylo=phylo,treedata=as.treedata(phylo),records=record)
          object@TREE<-TREEobject
          return(object)
          })


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create_mitoTracing
#'
#' This function is to create mitoTracing with basic information
#' @param GTsummary_list simply put GTSummary (Generated by CW_mgatk.read) into list, this allows mergeing multiple dataset this way.
#' @param depth_list  simply put depth(Generated by DepthSummary) into list, this allows mergeing multiple dataset this way.
#' @param feature.list_list  simply put feature.list(Generated by Vfilter_v3) into list, this allows mergeing multiple dataset this way.
#' @param labels a vector of labels for the samples.
#' @param thr One of the following "Total","VerySensitive","Sensitive","Specific"
#' @param qualifiedCellCut The minimum median mitochondrial coverage for a qualified cell, default is 10
#' @param OnlyHetero If only consider the heteroplasmy variants, default is T
#' @param VAFcut only use variants with VAF smaller than VAFcut. Default is 1.  We can use smaller value to constrain into only using rare variants
#' @param Cellcut only use variants with at least cellcut cells carry
#' @param maxctscut only use variants with at least in one cell with at leaset maxctscut variant fragments
#' @return mitoTracing class
#' @export
#' @import Seurat ape phytools phangorn treeio ggtree tidytree ggtreeExtra
Create_mitoTracing<-function(GTsummary_list,depth_list,feature.list_list,meta_list,labels,thr="VerySensitive",qualifiedCellCut=10,OnlyHetero=T,VAFcut=1,Cellcut=2,maxctscut=2){
require(ape)
require(phytools)
require(phangorn)
require(treeio)
require(ggtree)
require(tidytree)
require(ggtreeExtra)
CellMeta.all<-c()
GTsummary.all<-c()
V.union<-c()
V.fitered.list<-list()
len<-length(GTsummary_list)
for(i in 1:len){
CellMeta<-subset(depth_list[[i]]$Total[[2]],meanCov>=qualifiedCellCut)
names(CellMeta)[1]<-"Cell"
CellMeta<-merge(CellMeta,meta_list[[i]],by.x="Cell",by.y="ATACName")
CellMeta$Cell<-paste(CellMeta$Cell,i,sep="_")
CellMeta$Label<-labels[i]
GTsummary<-GTsummary_list[[i]][[thr]]
GTsummary$Cell<-paste(GTsummary$Cell,i,sep="_")
V.filtered<-subset(feature.list_list[[i]][[thr]],VAF<=VAFcut & CellN>=Cellcut & maxcts>=maxctscut)
if(OnlyHetero){
    V.filtered<-subset(V.filtered,HomoTag=="Hetero")
}
CellMeta.all<-rbind(CellMeta.all,CellMeta)
GTsummary.all<-rbind(GTsummary.all,GTsummary)
V.union<-c(V.union,as.character(V.filtered$Variants))
V.fitered.list<-c(V.fitered.list,list(V.filtered))
}
V.union.unique<-unique(V.union)
names(V.fitered.list)<-labels
GTsummary.all.filtered<-subset(GTsummary.all,Variants %in% V.union.unique & Cell %in% CellMeta.all$Cell)
ob<-mitoTracing()
ob@GTsummary.filtered<-GTsummary.all.filtered
ob@CellMeta<-CellMeta.all
ob@V.fitered.list=V.fitered.list
ob@UniqueV<-V.union.unique
ob@para<-c(Threhold=thr,qualifiedCellCut=qualifiedCellCut,OnlyHetero=OnlyHetero,VAFcut=VAFcut,Cellcut=Cellcut,maxctscut=maxctscut)
return(ob)
}


#' Compute distances for binary distances
#' @param M the binary matrix, Each row is a cell, each column is a variant, generated by Make_matrix
#' @param method distance method, choose from Jaccard, Dice, 3WJaccard, Simpson, Kulczynski2, Ochiai, Hamming
#' @examples d.Jaccard<-BinaryDist(object@Cts.Mtx.bi,method="Jaccard")
#' @export
#' @return dist object
BinaryDist<-function(M,method="Jaccard"){
print("This function compute pairwise distance(row-row) for binary matrix, input sparse matrix(Each row is cell, each column is variant)")
print("Available method:")
print(c("Jaccard","Dice","3WJaccard","Simpson","Kulczynski2","Ochiai","Hamming"))
require(Matrix)
Total<-Matrix::rowSums(M) # Compute total variant number for each cell
a<-M %*% Matrix::t(M)   ## Compute the overlaped variants across any two cells
b<-Total-a  ## Compute the variant only for the give row but not for the given column
c<-Matrix::t(b)  ## Compute the variant only for the give column but not for the given row
if(method=="Jaccard"){
disimilarity<-1-a/(a+b+c)
distance<-as.dist(disimilarity)
}else if(method=="Dice"){
disimilarity<-1-2*a/(2*a+b+c)
distance<-as.dist(disimilarity)
}else if(method=="Simpson"){
bcmin<-pmin(b,c)
disimilarity<-1-a/(a+bcmin)
distance<-as.dist(disimilarity)
}else if(method=="Kulczynski2"){
pr1<-a/(a+b)
pr2<-a/(a+c)
disimilarity<-1-(pr1+pr2)/2
distance<-as.dist(disimilarity)
}else if (method=="Ochiai"){
pr1<-a/(a+b)
pr2<-a/(a+c)
disimilarity<-1-sqrt(pr1*pr2)
distance<-as.dist(disimilarity)
}else if(method=="Hamming"){
disimilarity<-b+c
distance<-as.dist(disimilarity)
}else if(method=="3WJaccard"){
disimilarity<-1-3*a/(3*a+b+c)
distance<-as.dist(disimilarity)
}
return(distance)
}
