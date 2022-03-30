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
    jaccard3W="dist",
    w_jaccard="dist",
    w_cosine="dist",
    LSIdist="dist"
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
#' @slot Ctx.Mtx.depth A sparse matrix cell-mitoVariants(total counts for each position), store the variant count
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
            Ctx.Mtx.depth="matrix",
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

#' Add_DepthMatrix
#' Optional, add a matrix with same dimension with the Cts.Mtx and Cts.Mtx.bi, which display the depths
#' @param object mitoTracin class
#' @param QualifiedTotalCts a big source data, usually at XXX/mitoV/final
#' @export
setGeneric(name="Add_DepthMatrix", def=function(object,QualifiedTotalCts,...) standardGeneric("Add_DepthMatrix"))


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
#' @param d "jaccard" or "Dice" or "jaccard3W" or  "w_jaccard"  "w_cosine"  "LSIdist"
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
          definition=function(object,binary=T,res=0.6,lsidim=2:50,rmvariants=c("Variants310TC","Variants3109TC","Variants5764CT")){
          require(Signac)
          if(binary){
              Cts.Mtx.bi<-as.matrix(object@Cts.Mtx.bi)
              Cts.Mtx.bi<-Cts.Mtx.bi[,!colnames(Cts.Mtx.bi) %in% rmvariants]
              Cts.Mtx.bi<-Cts.Mtx.bi[rowSums(Cts.Mtx.bi)>0,]
              Cell_Variant.seurat<-CreateSeuratObject(counts = t(as.matrix(Cts.Mtx.bi)), assay = "mitoV")
          }else{
              Cts.Mtx<-as.matrix(object@Cts.Mtx)
              Cts.Mtx<-Cts.Mtx[,!colnames(object@Cts.Mtx) %in% rmvariants]
              Cts.Mtx<-Cts.Mtx[rowSums(Cts.Mtx)>0,]
              Cell_Variant.seurat<-CreateSeuratObject(counts = t(as.matrix(Cts.Mtx)), assay = "mitoV")
          }
          VariableFeatures(Cell_Variant.seurat) <- row.names(Cell_Variant.seurat) #names(which(Matrix::rowSums(Cell_Variant.seurat) > 100))
          Cell_Variant.seurat <- RunTFIDF(Cell_Variant.seurat, n = 50)
          Cell_Variant.seurat<- FindTopFeatures(Cell_Variant.seurat, min.cutoff = 2)
          Cell_Variant.seurat <- RunSVD(Cell_Variant.seurat, n = 50)
          Cell_Variant.seurat <- RunUMAP(Cell_Variant.seurat, reduction = "lsi", dims = lsidim)
          Cell_Variant.seurat <- FindNeighbors(Cell_Variant.seurat,reduction ="lsi"  ,dims = lsidim)
          Cell_Variant.seurat <- FindClusters(Cell_Variant.seurat, resolution = res)
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
#' @param jaccard  default=T
#' @param dice    default=T
#' @param jaccard3w  default=T
#' @param w_jaccard   default=T
#' @param w_cosine default=T
#' @param weight A two column dataframe, "Variant"(The variant name should match cell-variant matrix column, e.g, Variants310TC), "weight" (numeric)
#' @param NN To replace NA, which means a variant shown in the object is not shown in the weight vector, with a number, default is 1 for jaccard system. 
#' @param LSIdist default=T
#' @param dim the dimensions to use to calculate LSI distance default is 2:50
#' @return mitoTracing class
#' @export
setMethod(f="AddDist",
          signature="mitoTracing",
          definition=function(object,jaccard=T,dice=T,jaccard3w=T,w_jaccard=T,w_cosine=T,weightDF=NULL,NN=1,LSIdist=T,dim=2:50){
          d.Jaccard<-NA
          d.Dice<-NA    
          d.3WJaccard<-NA
          d.w_jaccard<-NA
          d.w_cosine<-NA
          if(length(weightDF)!=0){
            weight<-data.frame(Variants=colnames(object@Cts.Mtx.bi)) %>% merge(.,weightDF,by="Variants",all.x = T,sort = F) %>% .$weight
          }  
          if(length(which(is.na(weight)))!=0){
            weight[is.na(weight)]<-NN
            print("Some variant i weight is not found in cell-variant matrix, use 1")
          }
          if(length(weight)!=ncol(object@Cts.Mtx.bi)){
             stop("The length of weight does not match the variant numbers in the martix")
          }
          print("Weight vector matches well with the Cell-Variant matrix, continue...")
          if(jaccard){
              d.Jaccard<-BinaryDist(object@Cts.Mtx.bi,method="Jaccard")
          }    
          if(dice){
               d.Dice<-BinaryDist(object@Cts.Mtx.bi,method="Dice")
          }
          if(jaccard3w){
              d.3WJaccard<-BinaryDist(object@Cts.Mtx.bi,method="3WJaccard")
          }
          if(w_jaccard){
              if(length(weightDF)==0){
                  stop("Please input the weight, otherwise turn off the w_jaccard")
                  
              }
              d.w_jaccard<-quick_w_jaccard(object@Cts.Mtx.bi,w=weight)
          }
          if(w_cosine){
              if(length(weightDF)==0){
                  stop("Please input the weight, otherwise turn off the w_cosine")
              }
              d.w_cosine<-quick_w_cosine(object@Cts.Mtx.bi,w=weight)
          }
          if(LSIdist){
              d.lsi<-dist(object@Seurat@reductions$lsi@cell.embeddings[,dim])
          }
          object@DistObjects<-new("DistObjects",jaccard=d.Jaccard, Dice=d.Dice,jaccard3W=d.3WJaccard,w_jaccard=d.w_jaccard,w_cosine=d.w_cosine,LSIdist=d.lsi)
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



#' Add_DepthMatrix
#' Optional, add a matrix with same dimension with the Cts.Mtx and Cts.Mtx.bi, which display the depths
#' @param object mitoTracin class
#' @param QualifiedTotalCts a big source data, usually at XXX/mitoV/final,  If needed, edit V1, the cell name, which may have additional postfix due to combine
#' @return mitoTracing class
#' @export
#' @import reshape2
setMethod(f="Add_DepthMatrix",
          signature="mitoTracing",
          definition=function(object,QualifiedTotalCts){
          require(reshape2)
          colnames(QualifiedTotalCts)<-c("Cell","Pos","Total","VerySensitive","Sensitive","Specific")
          Dic<-gsub("Variants","",colnames(object@Cts.Mtx.bi)) %>% substr(.,1,nchar(.)-2) %>% as.integer %>% data.frame(Variants=colnames(object@Cts.Mtx.bi),Pos=.)
          QualifiedTotalCts.subset<-subset(QualifiedTotalCts,Cell %in% row.names(object@Cts.Mtx.bi)) %>% merge(.,Dic,by="Pos") %>% .[,c("Cell","Variants",object@para["Threhold"])]
          DepthMatrix<-dcast(QualifiedTotalCts.subset,Cell~Variants) %>% tibble::column_to_rownames("Cell") %>% as.matrix
          if (all(dim(object@Cts.Mtx.bi)==dim(DepthMatrix))){
             object@Ctx.Mtx.depth<-DepthMatrix[row.names(object@Cts.Mtx.bi),colnames(object@Cts.Mtx.bi)]
          }else{
              print(dim(object@Cts.Mtx.bi))
              print(dim(DepthMatrix))
              print("Check the input QualifiedTotalCts, the dimension cannot match")
          }
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


#' Compute weighted jaccard distance
#' @param M the binary matrix, Each row is a cell, each column is a variant, generated by Make_matrix
#' @param w weight for each variant, a vector
#' @export
#' @return dist object
quick_w_jaccard<-function(M,w){ 
## I tested against for loop based method, no problem    
    total<-M %*% w
    a<-M %*% (Matrix::t(M)*w)
    b<-as.numeric(total) - a
    c<-Matrix::t(b)
    disimilarity<-1-a/(a+b+c)
    distance<-as.dist(disimilarity)
    return(distance)
}

#' Compute weighted cosine distance
#' @param M the binary matrix, Each row is a cell, each column is a variant, generated by Make_matrix
#' @param w weight for each variant, a vector
#' @export
#' @return dist object
quick_w_cosine<-function(M,w){ 
mag<-as.numeric(sqrt(M %*% w))
a<-M %*% (Matrix::t(M)*w)
denominator<-outer(mag,mag)
disimilarity<-1-a/denominator
distance<-as.dist(disimilarity)
return(distance)
}



#' Make_AnnTable, Make a big dataframe, each row is a cell, each column includes info such as clonal UMAP, Clonal ID, ATAC/RNA/WNN UMAP, PCA, gene expression of chosen gene, etc.  Require a MitoTracing object and a multiome wrapper that better matches the cells in the MitoTracing  
#' @param Mitotracing  eg. DN4_HSC_mitoTracing.Sensitive
#' @param Multiome   eg. Donor04_HSC_Multiome_wrapper,  Multiome_wrapper object that matches with the MitoTracing, a reclustering using Multi_Wrapper() is recommended
#' @param clonal_features   eg. c("nCount_mitoV","seurat_clusters"), The column names take from Mitotracing@Seurat@meta.data, importantly the clonal clusterings
#' @param clonal_features_rename   eg. c("nCount_mitoV","clone_clusters") Rename the clonal_features
#' @param CellMeta_features   eg. c("meanCov","nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC","CellType") The column names take from  Mitotracing@CellMeta, may useful cell features
#' @param CellMeta_features_rename   eg. c("Mito_meanCov","nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC","CellType") Rename the CellMeta
#' @param multiome_features   eg. c("seurat_clusters")  The column names take from Multiome@meta.data 
#' @param multiome_features_rename   eg. c("NewSeurat_cluster")   Rename the column names for multiome_features
#' @param RNAUMAP    default T
#' @param ATACUMAP   Default T
#' @param WNNUMAP   Default T
#' @param PCA   Default T
#' @param LSI   Default T
#' @param genes   Default ""  can be a vector of gene names, for example c("HLF","CD34")
#' @param peaks   Default ""  can be a vector of peaks names
#' @param Variants Default ""  can be a vector of variant names format is eg "Variants10020TC"
#' @param PostTrans_from   Default c(2,3)  # This is a tricky part eh nmerging files are involved, find the postfix from cellranger agg for different sample
#' @param PostTrans_to   Default c(2,1)
#' @export
#' @import dplyr EZsinglecell2
#' @return AnnTable
Make_AnnTable<-function(
    Mitotracing=DN4_HSC_mitoTracing.Sensitive,
    Multiome=Donor04_HSC_Multiome_wrapper,
    clonal_features=c("nCount_mitoV","seurat_clusters"),
    clonal_features_rename=c("nCount_mitoV","clone_clusters"),
    CellMeta_features=c("meanCov","nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC","CellType"),
    CellMeta_features_rename=c("Mito_meanCov","nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC","CellType"),
    multiome_features=c("seurat_clusters"),
    multiome_features_rename=c("NewSeurat_cluster"),
    RNAUMAP=T,
    ATACUMAP=T,
    WNNUMAP=T,
    PCA=F,
    LSI=F,
    Variants="",
    genes="",
    peaks="",
    PostTrans_from=c(2,3),  # This is a tricky part eh nmerging files are involved, find the postfix from cellranger agg for different sample
    PostTrans_to=c(2,1)
){
multiome_meta_tb<-Translate_RNA2ATAC(Multiome@meta.data, PostFix =T, from = PostTrans_from, to=PostTrans_to)
multiome_meta_tb$RNAname<-row.names(multiome_meta_tb)
row.names(multiome_meta_tb)<-multiome_meta_tb$ATACName
## First, check if the cells in MitoTracing object are well matched in Multiome_wrapper
Matching<-length(intersect(Mitotracing@CellMeta$Cell,multiome_meta_tb$ATACName))
print(paste(length(Mitotracing@CellMeta$Cell),"Cells in MitoTracing object,", length(multiome_meta_tb$ATACName), "Cells in Multiome_wrapper object ---", Matching, "Cells matched"))
if(Matching/length(Mitotracing@CellMeta$Cell)<0.5){
    stop("Less than 10% of cells in MitoTracing is not matchable by multiome_wrapper, please check. Hint, maybe the PostTrans_from and to is messed up?")
}

## Make AnnTable, cellname|clone Umap1| clone umap2 (required)
AnnTable<-Mitotracing@Seurat@reductions$umap@cell.embeddings %>% as.data.frame %>%tibble::rownames_to_column("Cell")
row.names(AnnTable)<-AnnTable$Cell
names(AnnTable)<-c("Cell","cloUMAP_1","cloUMAP_2")
## Make Clonal feature table -- the clonal clustering, etc (Reccomended)
clonal_features_tb<-Mitotracing@Seurat@meta.data[,clonal_features,drop=F]
names(clonal_features_tb)<-clonal_features_rename
## Make CellMeta feature table -- mito coverage, ATAC/RNA counts, etc, these info are coming from initial multiome wrapping
CellMeta_features_tb<-Mitotracing@CellMeta[,CellMeta_features,drop=F]
names(CellMeta_features_tb)<-CellMeta_features_rename
## Make multiome feature table -- only to be useful when re-clustering is performed  
Multiome_feature_tb<-multiome_meta_tb[,multiome_features,drop=F]
names(Multiome_feature_tb)<-multiome_features_rename
## Any extra info can be put in, should be flexible. Such as signature(Cell cycle, etc) would be useful, which can be put in Multiome@meta.data
AnnTable<-Tomerge_v2(AnnTable,clonal_features_tb) %>% Tomerge_v2(.,CellMeta_features_tb) %>% Tomerge_v2(.,Multiome_feature_tb)

if(ATACUMAP){
umap.atac<-Translate_RNA2ATAC(as.data.frame(Multiome@reductions$umap.atac@cell.embeddings),PostFix =T, from = PostTrans_from, to=PostTrans_to) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("ATACName")
AnnTable<-Tomerge_v2(AnnTable,umap.atac)
}
if(RNAUMAP){
umap.rna<-Translate_RNA2ATAC(as.data.frame(Multiome@reductions$umap.rna@cell.embeddings),PostFix =T, from = PostTrans_from, to=PostTrans_to) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("ATACName")
AnnTable<-Tomerge_v2(AnnTable,umap.rna)
}
if(WNNUMAP){
wnn.umap<-Translate_RNA2ATAC(as.data.frame(Multiome@reductions$wnn.umap@cell.embeddings),PostFix =T, from = PostTrans_from, to=PostTrans_to) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("ATACName")
AnnTable<-Tomerge_v2(AnnTable,wnn.umap)
}
if(PCA){
pca<-Translate_RNA2ATAC(as.data.frame(Multiome@reductions$pca@cell.embeddings),PostFix =T, from = PostTrans_from, to=PostTrans_to) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("ATACName")
AnnTable<-Tomerge_v2(AnnTable,pca)
}
if(LSI){
lsi<-Translate_RNA2ATAC(as.data.frame(Multiome@reductions$lsi@cell.embeddings),PostFix =T, from = PostTrans_from, to=PostTrans_to) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("ATACName")
AnnTable<-Tomerge_v2(AnnTable,lsi)
}

if(length(genes)>0){
Available.genes<-genes[genes %in% row.names(Multiome@assays$SCT@data)]
expression<-Translate_RNA2ATAC(as.data.frame(t(as.matrix(Multiome@assays$SCT@data[Available.genes,]))),PostFix =T, from = PostTrans_from, to=PostTrans_to) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("ATACName")
AnnTable<-Tomerge_v2(AnnTable,expression)
}

if(length(peaks)>0){
Available.peaks<-peaks[peaks %in% row.names(Multiome@assays$ATAC@data)]
accessibility<-Translate_RNA2ATAC(as.data.frame(t(as.matrix(Multiome@assays$ATAC@data[Available.peaks,]))),PostFix =T, from = PostTrans_from, to=PostTrans_to) %>% tibble::remove_rownames() %>% tibble::column_to_rownames("ATACName")
names(accessibility)<-gsub("-","_",names(accessibility))
AnnTable<-Tomerge_v2(AnnTable,accessibility)
}
if(length(Variants)>0){
Available.Variants<-Variants[Variants %in% colnames(Mitotracing@Cts.Mtx.bi)]
variantstoshow<-as.matrix(Mitotracing@Cts.Mtx.bi[,Available.Variants])
AnnTable<-Tomerge_v2(AnnTable,variantstoshow)
}
AnnTable[is.na(AnnTable)]<-0
return(AnnTable)
}


#' Subset_MitoTracing Subset a mitotracing object by selecting a subset of cells, return a new MitoTracing object with only 4 slots: para; CellMeta; Cts.Mtx.bi; UniqueV, can be used for downstreme compute distance, clonal clustering, make tree, etc
#' @param Mitotracing  The Parent MitoTracing object eg. DN4_HSC_mitoTracing.Sensitive
#' @param Cells   Important, give a vector of Cell names(ATAC cell names)
#' @param ExtraInfo   Extra information, usually "Subset from ..."
#' @import Matrix
#' @export
#' @return MitoTracing Object
Subset_MitoTracing<-function(MitoTracing,Cells,ExtraInfo="Subset from ... "){
print("Raw slots are skipped: -GTsummary.filtered, -V.fitered.list, ")
Cells<-Cells[Cells %in% row.names(MitoTracing@Cts.Mtx.bi)]    
Variants.matrix.bi<-MitoTracing@Cts.Mtx.bi[Cells,]
Variants.matrix.bi.filtered<-Variants.matrix.bi[,colSums(Variants.matrix.bi)>=2]
Variants.matrix.bi.filtered<-Variants.matrix.bi.filtered[rowSums(Variants.matrix.bi.filtered)>=2,]  ## This step will filter out some cells
ob<-new("mitoTracing")
ob@para<-c(MitoTracing@para,Extra=ExtraInfo)
ob@CellMeta<-subset(MitoTracing@CellMeta,Cell %in% Cells)
ob@Cts.Mtx.bi<-Variants.matrix.bi.filtered
ob@UniqueV<-colnames(Variants.matrix.bi.filtered)
return(ob)
}

#' Define a function to make a list, each contains the cell names for a node 
#' @param tr phylo object (ape)
#' @param min.node.size default is 10, only the nodes with more than 10 tips are included  ( # Minimum # tips in the node to be included)
#' @param max.node.fra default is 0.33, only consider the nodes with less than max.node.fra*total cell number (# The up limit of the node size(Fraction of all tips) to be considered)
#' @import phangorn
#' @export
#' @return return a list each contains the cell names for a node that meets the criteria
Make_Cells4Nodes<-function(tr=DN4_SLCT_HSC_w_jaccard.njtree@phylo, min.node.size=10, max.node.fra=0.33){ 
require(phangorn)
max.node.size<-as.integer(length(tr$tip.label)*max.node.fra)
AllDes.list<-Descendants(tr, type = "tips")     
NodeSizes<-sapply(AllDes.list,length)
InterestNode<-which(NodeSizes> min.node.size& NodeSizes<max.node.size)
head(InterestNode)
Tips4Nodes.List<-Descendants(tr, node=InterestNode,type = "tips")
Cells4Nodes<-lapply(Tips4Nodes.List,function(x){tr$tip.label[x]})
names(Cells4Nodes)<-InterestNode    
return(Cells4Nodes)    
}