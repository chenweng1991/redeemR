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


#' Major redeem class that store clonal-resolved multi-omics
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
redeemR<-setClass(
    "redeemR",
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
            TREE="TREE",
            AssignedVariant="list"
           )
)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Make_matrix
#' This will make the matixies of Cell VS mitochondrial variants and return redeemR
#' Results stored in Cts.Mtx and Cts.Mtx.bi
#' @param object redeemR class
#' @export
setGeneric(name="Make_matrix", def=function(object)standardGeneric("Make_matrix"))

#' SeuratLSIClustering
#' This will use the mito variants for Seurat clustering (LSI based)
#' @param object redeemR class
#' @export
setGeneric(name="SeuratLSIClustering", def=function(object,...)standardGeneric("SeuratLSIClustering"))

#' AddDatatoplot_clustering
#' This prepare the clonal clustering data to plot
#' @param object redeemR class
#' @export
setGeneric(name="AddDatatoplot_clustering", def=function(object,...)standardGeneric("AddDatatoplot_clustering"))

#' AddDist
#' This add Jaccard, Dice, Jaccard3W distance and stored in DistObjects
#' @param object redeemR class
#' @export
setGeneric(name="AddDist", def=function(object,...) standardGeneric("AddDist"))

#' Make_tree
#' This will generate a basic phylogenetic tree
#' @param object redeemR class
#' @param d "jaccard" or "Dice" or "jaccard3W"
#' @param algorithm the algorithm used to build the tree, choose from "nj" and "upgma"
#' @export
setGeneric(name="Make_tree", def=function(object,d="jaccard", algorithm="upgma",onlyreturntree=F,...) standardGeneric("Make_tree"))

#' Add_Tree
#' Optional, if a phylogentic tree object phylo is already available, can be directly added to the redeemR
#' @param object redeemR class
#' @param phylo phyogenetic tree object
#' @export
setGeneric(name="AddTree", def=function(object,phylo,...) standardGeneric("AddTree"))

#' Add_DepthMatrix
#' Optional, add a matrix with same dimension with the Cts.Mtx and Cts.Mtx.bi, which display the depths
#' @param object redeemR class
#' @param QualifiedTotalCts a big source data, usually at XXX/mitoV/final
#' @export
setGeneric(name="Add_DepthMatrix", def=function(object,QualifiedTotalCts,...) standardGeneric("Add_DepthMatrix"))

#' Add_AssignVariant
#' a function to assign variants to edges based on maximum likihood
#' @param object redeemR class
#' @param QualifiedTotalCts a big source data, usually at XXX/mitoV/final
#' @export
setGeneric(name="Add_AssignVariant", def=function(redeemR,n.cores,...) standardGeneric("Add_AssignVariant"))

#' Add_tree_cut
#' a function to cut tree using assigned variant as branch-length on edge
#' @param redeemR  Need to have had the tree built
#' @param MinCell The minimum number of cells in each clone, otherwise merge with sibling
#' @param N  branch length to cut the tree
#' @export
setGeneric(name="Add_tree_cut", def=function(redeemR,MinCell,N,...) standardGeneric("Add_tree_cut"))



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Method definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' show
#' This will show the basics of redeemR class
#' @param object redeemR class
#' @return print out basics
#' @export
setMethod(f="show",
          signature="redeemR",
          definition=function(object){
              print(object@para)
              print(paste("Total Cell number:",nrow(object@CellMeta)))
              print(table(object@CellMeta$Label))
              print(paste("Total Variant number:",length(object@UniqueV)))
              print(paste("Slot:",slotNames(object)))
          })

#' Make_matrix
#' This will make the matixies of Cell VS mitochondrial variants and return redeemR
#' Results stored in Cts.Mtx and Cts.Mtx.bi
#' @param object mitoTracin class
#' @return mitoTracin class
#' @export
setMethod(f="Make_matrix",
          signature="redeemR",
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
          signature="redeemR",
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
#' @param  redeemR class
#' @param binary  Default is tree, to make use of the binary matrix
#' @param res     Default os 0.3, the resolution of the clustering
#' @return redeemR class
#' @export
setMethod(f="SeuratLSIClustering",
          signature="redeemR",
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
          Cell_Variant.seurat <- RunTSNE(Cell_Variant.seurat, reduction = "lsi", dims = lsidim,check_duplicates = FALSE)
          Cell_Variant.seurat <- FindNeighbors(Cell_Variant.seurat,reduction ="lsi"  ,dims = lsidim)
          Cell_Variant.seurat <- FindClusters(Cell_Variant.seurat, resolution = res)
          object@Seurat<-Cell_Variant.seurat
          return(object)
})

#' AddDatatoplot_clustering
#' This prepare the clonal clustering data to plot
#' @param object mitoTracin class
#' @return redeemR class
#' @export
setMethod(f="AddDatatoplot_clustering",
          signature="redeemR",
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
#' @return redeemR class
#' @export
setMethod(f="AddDist",
          signature="redeemR",
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
#' Optional, if a phylogentic tree object phylo is already available, can be directly added to the redeemR class in slot TREE
#' @param object mitoTracin class
#' @param phylo phyogenetic tree object
#' @param object mitoTracin class
#' @return redeemR class
#' @export
setMethod(f="AddTree",
          signature="redeemR",
          definition=function(object,phylo,record=""){
          TREEobject<-new("TREE",phylo=phylo,treedata=as.treedata(phylo),records=record)
          object@TREE<-TREEobject
          return(object)
          })



#' Add_DepthMatrix
#' Optional, add a matrix with same dimension with the Cts.Mtx and Cts.Mtx.bi, which display the depths
#' @param object mitoTracin class
#' @param QualifiedTotalCts a big source data, usually at XXX/mitoV/final,  If needed, edit V1, the cell name, which may have additional postfix due to combine
#' @return redeemR class
#' @export
#' @import reshape2
setMethod(f="Add_DepthMatrix",
          signature="redeemR",
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



#' a function to assign variants to edges based on maximum likihood
#' 
#' @param redeemR  Need to have redeemR@Ctx.Mtx.depth (By Add_DepthMatrix),  redeemR@Cts.Mtx  redeemR@Cts.Mtx.bi, redeemR@TREE
#' @export
#' @return redeemR with @AssignedVarian list of two p is a probability matrix of variants vs edges (Rowsum is 1) and Variant.assign.report, a dataframe (Variant|Edge.Assign|prob) 
#' @import foreach doParallel pryr doMC
setMethod(f="Add_AssignVariant",
          signature="redeemR",
          definition=function(redeemR=DN1_HSC_redeemR.VerySensitive,n.cores=4){
require(foreach)
# require(doParallel)
require(doMC)
require(pryr)
tree<-redeemR@TREE@phylo
mtr<-redeemR@Cts.Mtx %>%  as.matrix() %>% t  # Each row is a variant, each column is a cell
mtr.bi<-redeemR@Cts.Mtx.bi %>%  as.matrix() %>% t  # Each row is a variant, each column is a cell
depth<-redeemR@Ctx.Mtx.depth  %>% as.matrix() %>% t  # Each row is a variant, each column is a cell
mtr.bi.t<-t(mtr.bi) # each row is a cell
mtr.bi<-mtr.bi[,tree$tip.label]
## Prepare df_profile_mtx.t 
df<-reconstruct_genotype_summary(tree)
df_profile_mtx<-df2ProfileMtx(df)
df_profile_mtx<-matrix(as.numeric(df_profile_mtx),nrow=nrow(df_profile_mtx),ncol=ncol(df_profile_mtx)) # convert to numeric
df_profile_mtx.t<-t(df_profile_mtx)  # Each row is a cell, each column is a node
## Set expected global heteroplasmy level 
heteroplasmy<-mtr/depth 
Global.heteroplasmy<-heteroplasmy[!is.na(heteroplasmy)] %>% .[.!=0] %>% median  # this is a number ~ 0.03-0.05, an expected p for modeling dropout
## Generate a binomial-based probability distribution for each cell to model the dropout
Zero.p<-dbinom(0,depth,Global.heteroplasmy)
Zero.p[Zero.p==1]<-0.99
Zero.logV<-log(Zero.p)  # A matrix of log(probability of zero given w/ mutation in the cell)  Variant vs cell
NonZero.logV<-log(1-Zero.p) # A matrix of log(probability of 1 or >=1 given w/ mutation in the cell)  Variant vs cell
Zero.logP<-log(0.95) # A matrix of log(probability of zero given no mutation in the cell)  Variant vs cell
NonZero.logP<-log(1-0.95) # A matrix of log(probability of 1 or >=1 given no mutation in the cell)  Variant vs cell
print("(use doMC)Will gc in each loop; Befrore going into the loop, the memory use is:")
print(mem_used())
## Compute the loglik
# my.cluster <- parallel::makeCluster(n.cores,type="FORK")
# print(my.cluster)
# doParallel::registerDoParallel(cl = my.cluster)
registerDoMC(n.cores)
Loglik<-foreach (i =1:nrow(mtr.bi),.combine='rbind') %dopar%  # To loop through all variants
{
    ## Make indicator matries for 1-1 (Inside the clade, and is mutation) and 0-0 (outside of clade, and is not mutation) 
    x<-df_profile_mtx.t+mtr.bi[i,]
    x_11<-matrix(nrow=nrow(x),ncol=ncol(x))
    x_00<-matrix(nrow=nrow(x),ncol=ncol(x))
    x_11[which(x!=2)]<-0
    x_11[which(x==2)]<-1
    x_00[which(x!=0)]<- -1
    x_00[which(x==0)]<-1
    x_00[which(x_00==-1)]<-0
    ## Make indicator matries for 1-0 (Inside the clade, and is NOT mutation) and 0-1 (outside of clade, and is mutation) 
    x_10<-df_profile_mtx.t-mtr.bi[i,] 
    x_10[which(x_10==-1)]<-0
    x_01<-mtr.bi[i,]-df_profile_mtx.t
    x_01[which(x_01==-1)]<-0
    ## Compute the Loglikihood across all nodes
    Loglik_11<-x_11*NonZero.logV[i,]
    Loglik_10<-x_10*Zero.logV[i,]
    Loglik_01<-x_01*NonZero.logP
    Loglik_00<-x_00*Zero.logP
    Loglik.i<-colSums(Loglik_11+Loglik_10+Loglik_01+Loglik_00)
    ## Make the LogLik matrix Variant(row) verses nodes(columns)
    ## cleanup
    # rm(x,x_11,x_00,x_01,x_10,Loglik_10,Loglik_01,Loglik_00,Loglik_11)
    # gc()
    # # sink("Memory.record",append=TRUE)
    # cat(paste("i=",i,"\n",sep=""))
    # cat(paste(round(as.numeric(mem_used())/1000000000,2),"GB"))
    return(Loglik.i)
}
print("The foreach loop is completed, the memory use is")
mem_used()
row.names(Loglik)<-row.names(mtr.bi)
# parallel::stopCluster(cl = my.cluster)
## Compute the probability
n=dim(mtr.bi)[1]
edge_ml=apply(Loglik,1,which.max) %>% as.numeric
p=exp(Loglik-Loglik[(as.numeric(edge_ml)-1)*n+1:n])
p=p/rowSums(p)
# Make the report
Variant.assign.report<-data.frame(Edge.Assign=tree$edge[,1][apply(p,1,which.max)],prob=apply(p,1,max)) %>% .[order(.$prob,decreasing=T),]
redeemR@AssignedVariant<-list(p=p,Variant.assign.report=Variant.assign.report)
return(redeemR)
})


#' a function to cut tree using assigned variant as branch-length on edge
#' @param redeemR  Need to have had the tree built
#' @param MinCell The minimum number of cells in each clone, otherwise merge with sibling
#' @param N  branch length to cut the tree
#' @param Dumpcut Number of can be tolerated to be removed to fulfill the right side. The small value-> Less unassignment, big clones
#' @export
#' @return 
#' @import phangorn ape
setMethod(f="Add_tree_cut",
          signature="redeemR",
          definition=function(redeemR=DN4_stemcell_redeemR.seed.verysensitive,MinCell=30,N=1,prob.cut=0.3,Dumpcut=100){
require(phangorn)
phy=redeemR@TREE@phylo   
Allnodes=MakeAllNodes(redeemR,prob.cut=prob.cut)
AllAncestors<-Ancestors(phy,phy$edge[,2],type="all")
SumV.df_CloneNode<-data.frame(Node=phy$edge[,2],SumV=sapply(AllAncestors,function(x){sum(Allnodes[x,3])})) %>% subset(.,SumV>=N) ## Make a node table Node|SumV  
SumV.df_CloneNode.filtered<-SumV.df_CloneNode[Ancestors(phy,SumV.df_CloneNode$Node) %>% sapply(.,function(x){!any(x %in% SumV.df_CloneNode$Node)}),]  ## Leave only the relative most ancester
SumV.df_CloneNode.filtered<-data.frame(SumV.df_CloneNode.filtered,CladeSize=sapply(Descendants(phy,SumV.df_CloneNode.filtered$Node,type="tips"),length)) %>% .[order(.$CladeSize),]
FinalCloneNodes<-SumV.df_CloneNode.filtered[,c(1,3)] %>% .[order(.$CladeSize),] 
Smallclones<-c()
while(FinalCloneNodes$CladeSize[1]<MinCell){
    SibNode<-Siblings(phy,FinalCloneNodes$Node[1])
    if(length(Descendants(phy,SibNode)[[1]])<100){
        Parent<-data.frame(Node=Ancestors(phy,FinalCloneNodes$Node[1],type="parent"),CladeSize=sapply(Descendants(phy,Ancestors(phy,FinalCloneNodes$Node[1],type="parent")),length))
        FinalCloneNodes<-rbind(FinalCloneNodes[-1,],Parent)
        FinalCloneNodes<-FinalCloneNodes[order(FinalCloneNodes$CladeSize),]
    }else{
        Smallclones<-rbind(Smallclones,FinalCloneNodes[1,])
        FinalCloneNodes<-FinalCloneNodes[-1,]            
    }
}
FinalCloneNodes<-rbind(FinalCloneNodes,Smallclones)
FinalCloneNodes<-FinalCloneNodes[!duplicated(FinalCloneNodes),]
FinalCloneNodes<-FinalCloneNodes[Ancestors(phy,FinalCloneNodes$Node) %>% sapply(.,function(x){!any(x %in% FinalCloneNodes$Node)}),]
FinalCloneNodes<-FinalCloneNodes[order(FinalCloneNodes$CladeSize,decreasing=T),]
Bigclones<-c()
while(FinalCloneNodes$CladeSize[1]>100){
node=FinalCloneNodes$Node[1]
X<-as.list(rep(0,nrow(FinalCloneNodes)))
names(X)<-FinalCloneNodes$Node
CurrentSize<-sapply(Descendants(phy,node),length)
if(length(Descendants(phy,node,type="children"))>2){ 
    x<-Descendants(phy,node,type="children")
    x.size<-sapply(Descendants(phy,x),length)
    a<-Descendants(phy,node,type="children")[order(x.size,decreasing=T)[1]]
    b<-Descendants(phy,node,type="children")[order(x.size,decreasing=T)[2]]
}else{
    a<-Descendants(phy,node,type="children")[1]
    b<-Descendants(phy,node,type="children")[2]
}
a.size<-sapply(Descendants(phy,a),length)
b.size<-sapply(Descendants(phy,b),length)
if(min(a.size,b.size)>MinCell){
    subreturn_small<-data.frame(Node=c(a,b)[which.min(c(a.size,b.size))],CladeSize=sapply(Descendants(phy,c(a,b)[which.min(c(a.size,b.size))]),length))
    subreturn_big<-data.frame(Node=c(a,b)[which.max(c(a.size,b.size))],CladeSize=sapply(Descendants(phy,c(a,b)[which.max(c(a.size,b.size))]),length))
    FinalCloneNodes<-FinalCloneNodes[-1,]
    FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_small)
    FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_big)
}else{
    node_sub=c(a,b)[which.max(c(a.size,b.size))]
    a_sub<-Descendants(phy,node_sub,type="children")[1]
    b_sub<-Descendants(phy,node_sub,type="children")[2]
    a_sub.size<-sapply(Descendants(phy,a_sub),length)
    b_sub.size<-sapply(Descendants(phy,b_sub),length)
    if(min(a_sub.size,b_sub.size)>MinCell){
        subreturn_small.2<-data.frame(Node=c(a_sub,b_sub)[which.min(c(a_sub.size,b_sub.size))],CladeSize=sapply(Descendants(phy,c(a_sub,b_sub)[which.min(c(a_sub.size,b_sub.size))]),length))
        subreturn_big.2<-data.frame(Node=c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))],CladeSize=sapply(Descendants(phy,c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))]),length))
        FinalCloneNodes<-FinalCloneNodes[-1,]
        FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_small.2)
        FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_big.2)
        }
    else{
        Dump<-min(a_sub.size,b_sub.size)
        node_sub=c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))]
        while(Dump<Dumpcut){
          a_sub<-Descendants(phy,node_sub,type="children")[1]
          b_sub<-Descendants(phy,node_sub,type="children")[2]
          a_sub.size<-sapply(Descendants(phy,a_sub),length)
          b_sub.size<-sapply(Descendants(phy,b_sub),length)
          if(min(a_sub.size,b_sub.size)<MinCell & min(a_sub.size,b_sub.size)>0){
            Dump<-Dump+min(a_sub.size,b_sub.size)
            node_sub=c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))]
          }else{
            subreturn_small.2<-data.frame(Node=c(a_sub,b_sub)[which.min(c(a_sub.size,b_sub.size))],CladeSize=sapply(Descendants(phy,c(a_sub,b_sub)[which.min(c(a_sub.size,b_sub.size))]),length))
            subreturn_big.2<-data.frame(Node=c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))],CladeSize=sapply(Descendants(phy,c(a_sub,b_sub)[which.max(c(a_sub.size,b_sub.size))]),length))
            FinalCloneNodes<-FinalCloneNodes[-1,]
            FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_small.2)
            FinalCloneNodes<-rbind(FinalCloneNodes,subreturn_big.2)
            break
          }
        }
        if(Dump>=Dumpcut){
        FinalCloneNodes<-FinalCloneNodes[-1,]
        Bigclones<-rbind(Bigclones,data.frame(Node=node,CladeSize=sapply(Descendants(phy,node),length)))
        }
        }
 }
FinalCloneNodes<-FinalCloneNodes[order(FinalCloneNodes$CladeSize,decreasing=T),] 
}
FinalCloneNodes<-rbind(FinalCloneNodes,Bigclones)
FinalCloneNodes<-FinalCloneNodes[order(FinalCloneNodes$CladeSize,decreasing=T),] 
FinalCloneNodes<-FinalCloneNodes[complete.cases(FinalCloneNodes),]
FinalCloneNodes$Clone<-1:nrow(FinalCloneNodes)
SumV.df_CloneNode.filtered<-SumV.df_CloneNode.filtered[order(SumV.df_CloneNode.filtered$CladeSize,decreasing=T),]
SumV.df_CloneNode.filtered$Clone<-1:nrow(SumV.df_CloneNode.filtered)
FinalClone.report<-FinalCloneNodes %>% apply(.,1,function(x){data.frame(Cell=phy$tip.label[unlist(Descendants(phy,as.numeric(x[1])))],Clade_merge=x[1],Clone_merge=x[3],row.names = NULL)}) %>% do.call(rbind,.)
OriginClone.report<-SumV.df_CloneNode.filtered  %>% apply(.,1,function(x){data.frame(Cell=phy$tip.label[unlist(Descendants(phy,as.numeric(x[1])))],Clade=x[1],Clone=x[4],row.names = NULL)}) %>% do.call(rbind,.)
redeemR@CellMeta<-redeemR@CellMeta[,!colnames(redeemR@CellMeta) %in% c("Clade_merge","Clone_merge","Clade","Clone")]
redeemR@CellMeta<-merge(redeemR@CellMeta,merge(FinalClone.report,OriginClone.report,all=T),all.x=T)
return(redeemR)  # Note, may get NA, if only a few cells are NA, it is expected those are not assigned confidently
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Create_redeemR
#'
#' This function is to create redeemR with basic information
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
#' @return redeemR class
#' @export
#' @import Seurat ape phytools phangorn treeio ggtree tidytree ggtreeExtra
Create_redeemR<-function(GTsummary_list,depth_list,feature.list_list,meta_list,labels,thr="VerySensitive",qualifiedCellCut=10,OnlyHetero=T,VAFcut=1,Cellcut=2,maxctscut=2){
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
ob<-redeemR()
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



#' Make_AnnTable, Make a big dataframe, each row is a cell, each column includes info such as clonal UMAP, Clonal ID, ATAC/RNA/WNN UMAP, PCA, gene expression of chosen gene, etc.  Require a redeemR object and a multiome wrapper that better matches the cells in the redeemR  
#' @param redeemR  eg. DN4_HSC_redeemR.Sensitive
#' @param Multiome   eg. Donor04_HSC_Multiome_wrapper,  Multiome_wrapper object that matches with the redeemR, a reclustering using Multi_Wrapper() is recommended
#' @param clonal_features   eg. c("nCount_mitoV","seurat_clusters"), The column names take from redeemR@Seurat@meta.data, importantly the clonal clusterings
#' @param clonal_features_rename   eg. c("nCount_mitoV","clone_clusters") Rename the clonal_features
#' @param CellMeta_features   eg. c("meanCov","nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC","CellType") The column names take from  redeemR@CellMeta, may useful cell features
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
#' @import dplyr
#' @return AnnTable
Make_AnnTable<-function(
    redeemR=DN4_HSC_redeemR.Sensitive,
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
## First, check if the cells in redeemR object are well matched in Multiome_wrapper
Matching<-length(intersect(redeemR@CellMeta$Cell,multiome_meta_tb$ATACName))
print(paste(length(redeemR@CellMeta$Cell),"Cells in redeemR object,", length(multiome_meta_tb$ATACName), "Cells in Multiome_wrapper object ---", Matching, "Cells matched"))
if(Matching/length(redeemR@CellMeta$Cell)<0.5){
    stop("Less than 10% of cells in redeemR is not matchable by multiome_wrapper, please check. Hint, maybe the PostTrans_from and to is messed up?")
}

## Make AnnTable, cellname|clone Umap1| clone umap2 (required)
AnnTable<-redeemR@Seurat@reductions$umap@cell.embeddings %>% as.data.frame %>%tibble::rownames_to_column("Cell")
row.names(AnnTable)<-AnnTable$Cell
names(AnnTable)<-c("Cell","cloUMAP_1","cloUMAP_2")
## Make Clonal feature table -- the clonal clustering, etc (Reccomended)
clonal_features_tb<-redeemR@Seurat@meta.data[,clonal_features,drop=F]
names(clonal_features_tb)<-clonal_features_rename
## Make CellMeta feature table -- mito coverage, ATAC/RNA counts, etc, these info are coming from initial multiome wrapping
CellMeta_features_tb<-redeemR@CellMeta[,CellMeta_features,drop=F]
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
Available.Variants<-Variants[Variants %in% colnames(redeemR@Cts.Mtx.bi)]
variantstoshow<-as.matrix(redeemR@Cts.Mtx.bi[,Available.Variants])
AnnTable<-Tomerge_v2(AnnTable,variantstoshow)
}
AnnTable[is.na(AnnTable)]<-0
return(AnnTable)
}


#' Subset_redeemR Subset a redeemR object by selecting a subset of cells, return a new redeemR object with only 4 slots: para; CellMeta; Cts.Mtx.bi; UniqueV, can be used for downstreme compute distance, clonal clustering, make tree, etc
#' @param redeemR  The Parent redeemR object eg. DN4_HSC_redeemR.Sensitive
#' @param Cells   Important, give a vector of Cell names(ATAC cell names)
#' @param ExtraInfo   Extra information, usually "Subset from ..."
#' @import Matrix
#' @export
#' @return redeemR Object
Subset_redeemR<-function(redeemR,Cells,ExtraInfo="Subset from ... "){
print("Raw slots are skipped: -GTsummary.filtered, -V.fitered.list, ")
Cells<-Cells[Cells %in% row.names(redeemR@Cts.Mtx.bi)]    
Variants.matrix<-redeemR@Cts.Mtx[Cells,]
Variants.matrix.bi<-redeemR@Cts.Mtx.bi[Cells,]

Variants.matrix.filtered<-Variants.matrix[,colSums(Variants.matrix.bi)>=2]
Variants.matrix.bi.filtered<-Variants.matrix.bi[,colSums(Variants.matrix.bi)>=2]

Variants.matrix.bi.filtered<-Variants.matrix.bi.filtered[rowSums(Variants.matrix.bi.filtered)>=2,]  ## This step will filter out some cells
Variants.matrix.filtered<-Variants.matrix.filtered[rowSums(Variants.matrix.bi.filtered)>=2,]  ## This step will filter out some cells
Variants.matrix.filtered<-Variants.matrix.filtered[row.names(Variants.matrix.filtered) %in% row.names(Variants.matrix.bi.filtered),]
ob<-new("redeemR")
ob@para<-c(redeemR@para,Extra=ExtraInfo)
ob@CellMeta<-subset(redeemR@CellMeta,Cell %in% Cells)
Variants.matrix.filtered<-Variants.matrix.filtered[row.names(Variants.matrix.bi.filtered),]
ob@Cts.Mtx<-Variants.matrix.filtered
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

#' Define a function to make nn list, which can be further used to make adjacency matrix 
#' This scan row by row, looking for k.param nearest neighbours
#' @param d Distance matrix, can be a dist object or matrix 
#' @param k.param Default is 15
#' @export
#' @return return an nn list, which has two components: nn$idx and nn$dist
MakeNN<-function(d,k.param=15){
d<-as.matrix(d)
n.cells<-dim(d)[1]
object<-d
knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
knd.mat <- knn.mat 
for (i in 1:n.cells) { 
  knn.mat[i, ] <- order(object[i, ])[1:k.param] 
  knd.mat[i, ] <- object[i, knn.mat[i, ]] 
} 
nn.idx <- knn.mat[, 1:k.param]
nn.dist<- knd.mat[, 1:k.param]
nn<-list(idx=nn.idx,dist=nn.dist)
return(nn)
}


#' Define a function convert nn list to adjacency matrix that can be further used for igraph 
#' 
#' @param nn  nn list, which has two components: nn$idx and nn$dist
#' @export
#' @return return an nn.matrix. This is adjacency matrix can be input to igraph graph<-graph_from_adjacency_matrix(nn.matrix,diag = F,mode = "undirected")
NN2M<-function(nn){
k.param<-ncol(nn$idx) 
j <- as.numeric(x = t(x = nn$idx))
i <- ((1:length(x = j)) - 1) %/% k.param + 1
nn.matrix <- sparseMatrix(i = i, j = j, x = 1, dims = c(nrow(nn$idx), nrow(nn$idx)))
rownames(x = nn.matrix) <- rownames(x = nn$idx)
colnames(x = nn.matrix) <- rownames(x = nn$idx)
return(nn.matrix)
}



##----------------------2022-5-5 adpat treemut methods in to assign variants into edge
#' This is a function borrowed from https://github.com/NickWilliamsSanger/treemut/blob/main/R/treemut.R#L68
#' Input phylo object, return a "profile matrix"--Edge(or denoted as the ending node) vs cell. a 0, 1 character string that indicate what cells in a given node
#' 
#' @param phylo  phylo an ape object
#' @export
#' @return df  includes df$df which is a big data frame,  and df$sample that is the cell names
reconstruct_genotype_summary=function(phylo
){
  dat=phylo$edge
  samples=phylo$tip.label
  N=length(samples)
  zeros=rep(0,N)
  profile=sapply(1:length(samples),function(i){tmp=zeros;tmp[i]=1;paste(tmp,collapse="")})
  df=data.frame(profile=profile,edge_length=phylo$edge.length[1:N])
  #Create empty list where each element will correspond to a node
  muts=lapply(1:dim(dat)[1],function(i) c())
  map_node=match(1:max(dat[,2]),dat[,2])
  #Here we loop through the tips (samples) and get all the ancestral nodes for each tip.
  #Then for each node in the list we add the tip to that nodes list of tips, ultimately giving us a list of samples that share the node. 
  for(i in 1:N){
    parents=get_ancestral_nodes(i,edge = dat,exclude_root=TRUE)
    for(j in parents){
      muts[[map_node[j]]]=append(muts[[map_node[j]]],i)
    }
  }
  profiles=sapply(muts,function(x){tmp=zeros;tmp[x]=1;paste(tmp,collapse="")})
  df=data.frame(profile=profiles,edge_length=phylo$edge.length,stringsAsFactors=FALSE)
  df=add_derived_profile_info(df,phylo$tip.label)
  list(df=df,samples=phylo$tip.label)
}
              
#' This is a convinience function, internal 
str2vector<-function(x){
    strsplit(x,"")
}
#' This is a convinience function, internal
df2ProfileMtx<-function(df){
    profile_mtx<-str2vector(df$df[,1,drop=T]) %>% do.call(rbind,.)
    return(profile_mtx)
}

#' This is a convinience function, internal borrowed from https://github.com/NickWilliamsSanger/treemut/blob/main/R/treemut.R#L68
get_ancestral_nodes=function(node,edge,exclude_root=TRUE){
  idx=which(edge[,2]==node)
  parents=node ##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=edge[idx,1]
    parents=c(parents,parent)
    #This finds the parent of the current parent - thus navigating up to the root.
    idx=which(edge[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)] ##The last node is the root.
  }else{
    parents
  }
}


#' This is a convinience function, internal borrowed from https://github.com/NickWilliamsSanger/treemut/blob/main/R/treemut.R#L68
add_derived_profile_info=function(profile_df,samples=sprintf("s%s",0:(nchar(profile_df$profile[1])-1))){
  base=rep(0,nchar(profile_df$profile[1]))
  samples_private=sapply(1:length(base),function(i){this=base;this[i]=1;paste(this,collapse="")})
  missing_samples=setdiff(samples_private,profile_df$profile)
  if(length(missing_samples)>0){
    profile_df=rbind(profile_df,data.frame(profile=missing_samples,edge_length=0))
  }
  
  profile_df$mut_count=nchar(profile_df$profile)-nchar(gsub("1","",profile_df$profile))
  profile_df$profile_int=lapply(profile_df$profile,split_profile)
  profile_df$id=1:dim(profile_df)[1]
  profile_df$label=profile_df$profile
  idx=which(profile_df$mut_count==1)
  profile_df$label[idx]=sapply(idx,function(i) samples[which(profile_df$profile_int[[i]]==1)])
  profile_df
}
#' This is a convinience function, internal  borrowed from https://github.com/NickWilliamsSanger/treemut/blob/main/R/treemut.R#L68
split_profile=function(profile){
  as.integer(unlist(strsplit(profile,split="")))
}


#' Define a function make the Allnodes(Node|Parent|Freq|CladeSize), where Freq is the number of variants assigned to the node(as ending point) 
#' from redeemR object, 
#' 
#' @param redeemR  a redeemR object already have the tree built
#' @param prob.cut The probability cutoff to include confidently assigned variant
#' @export
#' @import dplyr
MakeAllNodes<-function(redeemR=DN4_stemcell_redeemR.seed.verysensitive, prob.cut=0.3){
phy=redeemR@TREE@phylo
NodeVariant.summary<-redeemR@AssignedVariant$Variant.assign.report %>% subset(prob>=prob.cut) %>% .$Edge.Assign %>% table %>% as.data.frame() %>% .[order(.$Freq,decreasing=T),]
colnames(NodeVariant.summary)[1]<-"Node"
EdgeTable<-phy$edge %>% as.data.frame
colnames(EdgeTable)<-c("Parent","Node")
Allnodes<- merge(EdgeTable,NodeVariant.summary,by="Node",all.x=T)
Allnodes$CladeSize<-Descendants(phy,Allnodes$Node) %>% sapply(.,length)
Allnodes$Freq[is.na(Allnodes$Freq)]<-0
return(Allnodes)
}


#' Define a function to perform single-cell based hard porogeny assignment
#' This function was developed based on DN4T2.basics.ipynb
#' @param HSC_redeemR  The HSC_redeemR is the redeemR object for defined HSC
#' @param Full_redeemR The FULL_redeemR is the redeemR object for the full BMMC_HSPC_HSC
#' @param distCut Default is 0.95, the distance, below which I define as the related progeny
#' @export
#' @import ggnewscale tibble dplyr RColorBrewer
ProgenyMapping<-function(HSC_redeemR=DN4_PhenoHSC_redeemR.verysensitive,Full_redeemR=DN4_BMMC_HSPC_HSC_redeemR.verysensitive,distCut=0.95,d="w_jaccard"){
require(ggnewscale)
require(tibble)
require(dplyr)
require(RColorBrewer)
## Prepare distance matrix from HSC to non-HSCs
BMMC_HSPC_HSC.dist<-as.matrix(slot(Full_redeemR@DistObjects,d))
HSCs<-HSC_redeemR@CellMeta$Cell
HSCs<-HSCs[HSCs %in% colnames(BMMC_HSPC_HSC.dist)]
NonHSCs<-row.names(BMMC_HSPC_HSC.dist)[!row.names(BMMC_HSPC_HSC.dist) %in% HSCs]
BMMC_HSPC_HSC.dist.fromHSC<-BMMC_HSPC_HSC.dist[NonHSCs,HSCs]

p1.1<-BMMC_HSPC_HSC.dist.fromHSC[,sample(1:ncol(BMMC_HSPC_HSC.dist.fromHSC),1)] %>% data.frame(Distance=.,rk=rank(.)) %>% ggplot()+aes(rk,Distance)+geom_point()+geom_hline(yintercept = distCut)+theme_bw()
p1.2<-BMMC_HSPC_HSC.dist.fromHSC[,sample(1:ncol(BMMC_HSPC_HSC.dist.fromHSC),1)] %>% data.frame(Distance=.,rk=rank(.)) %>% ggplot()+aes(rk,Distance)+geom_point()+geom_hline(yintercept = distCut)+theme_bw()
p1.3<-BMMC_HSPC_HSC.dist.fromHSC[,sample(1:ncol(BMMC_HSPC_HSC.dist.fromHSC),1)] %>% data.frame(Distance=.,rk=rank(.)) %>% ggplot()+aes(rk,Distance)+geom_point()+geom_hline(yintercept = distCut)+theme_bw()
meta<-Full_redeemR@CellMeta[,c("Cell","meanCov")] %>% tibble::column_to_rownames("Cell")
HSC.output<-apply(BMMC_HSPC_HSC.dist.fromHSC,2,function(x){length(which(x<distCut))}) %>% data.frame(OutputN=.) %>% Tomerge_v2(.,meta)
VariantNumber<-Full_redeemR@GTsummary.filtered %>% subset(.,Cell %in% HSCs) %>% group_by(Cell) %>% dplyr::summarise(VariantNumber=n()) %>% tibble::column_to_rownames(var="Cell")
HSC.output<-Tomerge_v2(HSC.output,VariantNumber)
HSC.output<-HSC.output %>% tibble::rownames_to_column(var = "Cell")
p2<-ggplot(HSC.output)+aes(VariantNumber,log10(OutputN))+geom_point()+xlim(0,40)+theme_bw()+geom_smooth(method="lm")
## Make a dic to translate from cell name to coarse cell type
dic<-recode(Full_redeemR@CellMeta$STD.CellType,
HSC="Stem",
MPP="EarlyP",CMP="EarlyP",MKP="MK",
MEP="ME",GMP="Mye",MDP="Mye",CDP="Mye", EryP="ME",
LMPP="Lym",CLP="Lym",ProB="Lym",Plasma="Lym",
Mono="Mye",Ery="ME",mDC="Mye", cDC="Mye",
CD4="Lym",CD8="Lym" ,NK="Lym"  ,B="Lym" ,pDC="Lym") 
names(dic)<-Full_redeemR@CellMeta$Cell
## Calculate output lineages
LineageNames<-levels(dic)
ResolvedLineage<-function(x){
outputLineage<-dic[names(x[x<0.95 & x>0])] %>% table %>% as.numeric 
return(outputLineage)
}
HSC_output_lineages<-apply(BMMC_HSPC_HSC.dist.fromHSC,2,ResolvedLineage)  %>% t() %>% as.data.frame #%>% rownames_to_column(var="Cell")
colnames(HSC_output_lineages)<-c(LineageNames)
HSC_output_lineages$Totaloutput<-rowSums(HSC_output_lineages)
## Expected lineage outputs
Lym.expM<-sum(HSC_output_lineages$Lym)/sum(HSC_output_lineages$Totaloutput)
Mye.expM<-sum(HSC_output_lineages$Mye)/sum(HSC_output_lineages$Totaloutput)
MK.expM<-sum(HSC_output_lineages$MK)/sum(HSC_output_lineages$Totaloutput)
ME.expM<-sum(HSC_output_lineages$ME)/sum(HSC_output_lineages$Totaloutput)
EarlyP.expM<-sum(HSC_output_lineages$EarlyP)/sum(HSC_output_lineages$Totaloutput)
## Define a function to analyze lineage  bias for each HSC
LineageBiasTest<-function(df=HSC_output_lineages,L="Lynphoid",L.exp=Lymphoid.expM)
{
Fold<-(df[,L]/df$Totaloutput)/L.exp
ps<-c()
for(i in 1:nrow(df)){
if(df$Totaloutput[i]==0){
    ps<-c(ps,1)
}else{
mod<-binom.test(df[,L][i],df$Totaloutput[i],L.exp,alternative = "greater")
ps<-c(ps,mod$p.value)
}
}
df<-data.frame(row.names=row.names(df),Fold,ps) 
return(df)
}
## Compute lineage bias fold change and p value
LineageBiasTest.Lym<-LineageBiasTest(df=HSC_output_lineages,L="Lym",L.exp=Lym.expM) %>% tibble::rownames_to_column(var="Cell")
LineageBiasTest.Mye<-LineageBiasTest(df=HSC_output_lineages,L="Mye",L.exp=Mye.expM)  %>% tibble::rownames_to_column(var="Cell")
LineageBiasTest.MK<-LineageBiasTest(df=HSC_output_lineages,L="MK",L.exp=MK.expM)  %>% tibble::rownames_to_column(var="Cell")
LineageBiasTest.ME<-LineageBiasTest(df=HSC_output_lineages,L="ME",L.exp=ME.expM)  %>% tibble::rownames_to_column(var="Cell")
LineageBiasTest.EarlyP<-LineageBiasTest(df=HSC_output_lineages,L="EarlyP",L.exp=EarlyP.expM)  %>% tibble::rownames_to_column(var="Cell")
n=length(unique(HSC_redeemR@CellMeta$Clone_merge))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
meta<-HSC_redeemR@CellMeta
meta$Clone_merge<-as.factor(meta$Clone_merge)
options(repr.plot.width=10,repr.plot.height=4,repr.plot.res=600)
p<-ggtree(HSC_redeemR@TREE@treedata,size=0.1,branch.length = "none")#+geom_tiplab(as_ylab=TRUE, color='firebrick')
p3<-p+
geom_facet(panel="Clone",data=meta,geom = geom_bar,mapping = aes(x = 1,fill=Clone_merge),orientation = 'y', width = 1, stat='identity')+
scale_fill_manual(values=col_vector)+
new_scale_fill()+
geom_facet(panel="Outputs",data=HSC.output,geom = geom_bar,mapping = aes(x = OutputN),orientation = 'y', width = 1, stat='identity',fill="darkgreen")+
geom_facet(panel="coverage",data=HSC.output,geom = geom_bar,mapping = aes(x = meanCov),orientation = 'y', width = 1, stat='identity',fill="gold3")+
geom_facet(panel="Lym",data=LineageBiasTest.Lym,geom = geom_bar,mapping = aes(x = log2(Fold),fill=-log10(ps)),orientation = 'y', width = 1, stat='identity')+
geom_facet(panel="Mye",data=LineageBiasTest.Mye,geom = geom_bar,mapping = aes(x = log2(Fold),fill=-log10(ps)),orientation = 'y', width = 1, stat='identity')+
geom_facet(panel="MK",data=LineageBiasTest.MK,geom = geom_bar,mapping = aes(x = log2(Fold),fill=-log10(ps)),orientation = 'y', width = 1, stat='identity')+
geom_facet(panel="ME",data=LineageBiasTest.ME,geom = geom_bar,mapping = aes(x = log2(Fold),fill=-log10(ps)),orientation = 'y', width = 1, stat='identity')+
geom_facet(panel="EarlyP",data=LineageBiasTest.EarlyP,geom = geom_bar,mapping = aes(x = log2(Fold),fill=-log10(ps)),orientation = 'y', width = 1, stat='identity')+
scale_fill_gradient(low="grey",high="red",limits=c(0.5,3),oob=scales::squish)+theme_tree2() +xlim_expand(c(-2,2), "Lym")+xlim_expand(c(-2,2), "Mye")+
xlim_expand(c(-2,2), "MK")+xlim_expand(c(-2,2), "ME")+xlim_expand(c(-2,2), "EarlyP")+theme(axis.text=element_text(color="black"))#+facet_widths(p3, widths=c(2,0.25,1,1,1,1,1,1,1)) #+theme_tree2(legend.position=c(.05, .85))
p3<-facet_widths(p3, widths=c(2,0.25,1,1,1,1,1,1,1))
#To meake a lineage summary, 
#output_lineages<-HSC_T2_output_lineages
#HSC_redeemR ## It has to be a redeemR object with Clone_merge
output_lineage.summary<-HSC_output_lineages %>% tibble::rownames_to_column("Cell") %>% merge(.,HSC_redeemR@CellMeta[,c("Cell","Clone_merge")]) %>% 
group_by(Clone_merge) %>% dplyr::summarise(Lym=sum(Lym),Mye=sum(Mye),ME=sum(ME),MK=sum(MK),TotalOutPut=sum(Totaloutput)) %>% subset(.,!is.na(Clone_merge)) %>% 
merge(.,as.data.frame(table(HSC_redeemR@CellMeta$Clone_merge)),by.x="Clone_merge",by.y="Var1") %>% mutate(OutLevel=TotalOutPut/Freq)
output_lineage.summary.pct<-output_lineage.summary[,1:5] %>% tibble::column_to_rownames("Clone_merge") %>% apply(.,1,function(x){x/sum(x)}) %>% t() 
output_lineage.summary.pct.scale<-scale(output_lineage.summary.pct) %>% data.frame(.,OutLevel.scale=scale(output_lineage.summary$OutLevel)) %>%
tibble::rownames_to_column("Clone_merge")
## Make list of analyzed dataframe to output
LineageSummary<-list(output_lineage.summary=output_lineage.summary,output_lineage.summary.pct=output_lineage.summary.pct,output_lineage.summary.pct.scale=output_lineage.summary.pct.scale)
LineageEnrich<-list(LineageBiasTest.Lym=LineageBiasTest.Lym,LineageBiasTest.Mye=LineageBiasTest.Mye,LineageBiasTest.MK=LineageBiasTest.MK,LineageBiasTest.ME=LineageBiasTest.ME,LineageBiasTest.EarlyP=LineageBiasTest.EarlyP)
return(list(LineageSummary=LineageSummary,LineageEnrich=LineageEnrich,p1_list=list(p1.1,p1.2,p1.3),p2=p2,p3=p3))
}

#' Define a function to perform Find marker for top vs bottom clones
#' This function was developed based on DN4T2.basics.ipynb
#' @param topClones  a vector of clone ID eg. c("1","3","7"),this must be in HSC_redeemR@CellMeta$Clone_merge
#' @param bottomClones  a vector of clone ID eg. c("2","5"), this must be in HSC_redeemR@CellMeta$Clone_merge
#' @param ob Seurat object (Multiomics), the postfix needs to be compatible with HSC_redeemR, the cells will be matched by cell names
#' @param HSC_redeemR  redeemR object for HSC
#' @param test  the statistic method to use for DE, a wrapper function from Seurat FindAllMarkers
#' @export 
#' @import tibble Seurat
Clone_FinderMarker<-function(topClones,bottomClones,HSC_Multiome_wrapper=Donor04_HSC_Multiome_wrapper,HSC_redeemR,assay="SCT",test="wilcox"){
require(Seurat)
clone.df<-HSC_redeemR@CellMeta[,c("Cell","Clone_merge")] %>% mutate(Cell=Translate_simple_ATAC2RNA(Cell)) %>% tibble::column_to_rownames("Cell")
meta.top<-subset(clone.df,Clone_merge %in% topClones) %>% mutate(LB="top")
meta.bottom<-subset(clone.df,Clone_merge %in% bottomClones) %>% mutate(LB="bottom")
meta<-rbind(meta.top,meta.bottom)
HSC_Multiome_wrapper@meta.data<-Tomerge_v2(HSC_Multiome_wrapper@meta.data,meta)
HSC_Multiome_wrapper<- subset(x = HSC_Multiome_wrapper,subset= LB %in% c("top","bottom"))
HSC_Multiome_wrapper<- SetIdent(HSC_Multiome_wrapper, value = "LB")
DefaultAssay(HSC_Multiome_wrapper) <- assay
markers <- FindAllMarkers(HSC_Multiome_wrapper, only.pos=T,min.pct = 0.01, logfc.threshold = 0, test.use=test,return.thresh = 0.99) # ,test.use="DESeq2"
return(markers)
}


#' FromDist2Graph 
#' From disttance object or matrix to graph, default is to return igraph object
#' This function was developed based on 
#' @param d the distance matrix,  this can be either dist or a matrix
#' @param k.param K default is 30
#' @param return_igraph Wheather return igraph, default is T which return igraph. Otherwise, return adjacent matrix
#' @return igraph or adjacent matrix
#' @export
#' @import Matrix igraph
FromDist2Graph<-function(d,k.param=30,return_igraph=T){
if(!is.matrix(d)){
  d<-as.matrix(d)
}
n.cells<-dim(d)[1]
knn.mat <- matrix(data = 0, ncol = k.param, nrow = n.cells)
knd.mat <- knn.mat 
for (i in 1:n.cells) { 
  knn.mat[i, ] <- order(d[i, ])[1:k.param] 
  knd.mat[i, ] <- d[i, knn.mat[i, ]] 
} 
nn.dist <- knn.mat[, 1:k.param]
knn<-nn.dist
# nn.dist<- knd.mat[, 1:k.param]
# nn<-list(idx=nn.idx,dist=nn.dist)
# Compute MNN  -- Borrowed from SCAVENGE https://github.com/sankaranlab/SCAVENGE
knn2 <- list()
    length(knn2) <- nrow(knn)
    for (i in 1:nrow(knn)) {
        xx <- apply(knn[knn[i, -1], ], 1, function(x) {
            any(x == i)
        })
        if (sum(xx) > 0) {
            temp_knn <- knn[i, c(TRUE, xx)]
            temp_el <- cbind(temp_knn[1], c(temp_knn[-1]))
        }
        else {
            temp_el <- knn[i, 1:2]
        }
        knn2[[i]] <- temp_el
    }
el <- do.call(rbind.data.frame, knn2) %>% as.matrix
adj <- igraph::get.adjacency(igraph::graph.edgelist(el))
mutualknn <- 1 * ((adj + Matrix::t(adj)) > 0)
colnames(mutualknn) <- rownames(mutualknn) <- rownames(d)
if(return_igraph){
  g<-graph_from_adjacency_matrix(mutualknn,diag = F,mode = "undirected")
  return(g)
}else{
  return(mutualknn)
}
}





#' ProgenyMapping_np 
#' Define a function to compute network propagation based probability
#' FromDist2Graph is needed to convert fistance matrix into MNN graph
#' @param HSC_redeemR  The HSC_redeemR is the redeemR object for defined HSC,  have already gone through Add_DepthMatrix--Add_AssignVariant--Add_tree_cut, otherwise, need
#' othereise, need a column in CellMeta that indicates the clone ID
#' @param Full_redeemR The FULL_redeemR is the redeemR object for the full BMMC_HSPC_HSC
#' @param CloneCol "Clone_merge"
#' @param k  the k.param used for MNN graph 
#' @param gm gamma default is 0.05 which mean 95% information is passing out
#' @param ProbCut The cutoff of the maximum probability for a given progeny cell(If the maximum probability is lower than ProbCut, it will be filtered)
#' @param Celltype The column to be used in aggregate into lineages
#' @return a list of two  ALLmeta.npClone (A meta data with last column npClone), np_mat (the network propagation matrix))
#' @export
#' @import SCAVENGE dplyr
ProgenyMapping_np<-function(HSC_redeemR=DN4_stemcell_redeemR.seed.verysensitive,
                       Full_redeemR=DN4_BMMC_HSPC_HSC_redeemR.verysensitive,
                       CloneCol="Clone_merge",k=30,gm=0.5,useLSI=F,useSCAVENGE_LSI=F,subsample=F,ProbCut=0.7,Celltype="Rig.CellType"){
print("Note: By default, usLSI=F, the MNN nretwork is based on jaccard; Alternatively, useLSI=T")
print("Note: In case useLSI,useSCAVENGE_LSI=T to use SCAVENGE to compute LSI, else, do LSI via Seurat")
print("Note: subsample=T, to subsample into same number of seed cells in each clone")
require(SCAVENGE)
require(dplyr)
HSCmeta<-HSC_redeemR@CellMeta
HSCmeta<-HSCmeta[!is.na(HSCmeta[,CloneCol]),]
ALLmeta<-Full_redeemR@CellMeta
if(subsample){
  print(table(HSCmeta[,CloneCol]))
  HSCmeta<-HSCmeta %>% group_by(Clone_merge) %>% sample_n(min(as.numeric(table(HSCmeta[,CloneCol]))))
  print(table(HSCmeta[,CloneCol]))
}
if(useLSI){
  print("use LSI")
    if(useSCAVENGE_LSI){
      print("use LSI-SCAVENGE")  
      rm_v = c("Variants310TC","Variants3109TC","Variants5764CT")
      rm_idx <- which(colnames(Full_redeemR@Cts.Mtx.bi) %in% rm_v)
      cell_rm_idx <- which(rowSums(Full_redeemR@Cts.Mtx.bi[, -rm_idx]) ==0)
      tfidf_mat = tfidf(bmat=t(Full_redeemR@Cts.Mtx.bi[-cell_rm_idx, -rm_idx]), mat_binary=TRUE, TF=TRUE, log_TF=TRUE)
      lsi_mat <- do_lsi(tfidf_mat, dims=30)
      AdjMtx<- SCAVENGE::getmutualknn(lsi_mat[, 2:k], k)
    }else{
    print("use LSI-Seurat")  
    AdjMtx<- SCAVENGE::getmutualknn(Full_redeemR@Seurat@reductions$lsi@cell.embeddings[, 2:k], k)
    }
}else{
  print("use raw w_jaccard")
AdjMtx<-FromDist2Graph(Full_redeemR@DistObjects@w_jaccard,k.param = k,return_igraph = F)
}
ALLCells<-row.names(AdjMtx)
np_mat <- matrix(nrow=length(ALLCells), ncol=length(table(HSCmeta[,CloneCol])))
rownames(np_mat) <- ALLCells
colnames(np_mat) <- names(table(HSCmeta[,CloneCol]))

for (i in 1:ncol(np_mat)){
    print(i)
    idx_temp <- HSCmeta[,CloneCol]==names(table(HSCmeta[,CloneCol]))[i]
    np_mat[, i] <- randomWalk_sparse(intM=AdjMtx, HSCmeta$Cell[idx_temp], gamma=gm)

}
ALLmeta.npClone<-data.frame(Cell=row.names(np_mat),npClone=colnames(np_mat)[apply(np_mat,1,which.max)]) %>% merge(ALLmeta[,c("Cell","meanCov","nCount_RNA","nCount_ATAC","STD.CellType","Rig.CellType","STD_Cat","STD_Cat2","Label")],.)
## Updated 2022-6-5 evaluate the probability
Probability<-apply(np_mat,1,function(x){x/sum(x)}) %>% as.data.frame()
Probability[is.na(Probability)]<-0
MaxProbability<-apply(Probability,2,max) %>% data.frame(MaxProb=.)
row.names(ALLmeta.npClone)<-ALLmeta.npClone$Cell
ALLmeta.npClone<-Tomerge_v2(ALLmeta.npClone,MaxProbability)
ALLmeta.npClone.filtered<-subset(ALLmeta.npClone,MaxProb>=ProbCut)
## To evaluate the HSC assignment
HSC.assign<-subset(ALLmeta.npClone,STD.CellType=="HSC") %>% .[,c("Cell","npClone","MaxProb")]
HSC.ori<-HSC_redeemR@CellMeta[,c("Cell","Clone_merge")] %>% .[complete.cases(.),] 
HSC.ori_assign<-merge(HSC.ori,HSC.assign)
Accuracy_HSC.BenchMark<-length(which(HSC.ori_assign$Clone_merge==HSC.ori_assign$npClone))/nrow(HSC.ori_assign)


ClonalSize<-as.data.frame(table(HSC_redeemR@CellMeta$Clone_merge)) %>% rename(npClone="Var1")
datatoplot<-cbind(ALLmeta.npClone.filtered[,c("Cell","npClone")],Lin=recode(ALLmeta.npClone.filtered[,Celltype],
HSC="Stem",
MPP="EarlyP",CMP="EarlyP",
MKP="MK",
MEP="ME",GMP="Mye",MDP="Mye", EryP="ME",
LMPP="EarlyP",CLP="Lym",ProB="Lym",Plasma="Lym",
Mono="Mye",Ery="Mye",cDC="Mye",
CD4="Lym",CD8="Lym" ,NK="Lym" ,B="Lym" ,pDC="Lym",ProB="Lym")) %>% group_by(npClone) %>% reshape2::dcast(npClone~Lin,fun.aggregate = length) %>% mutate(Total=EarlyP+Lym+ME+MK+Mye) %>% 
merge(.,ClonalSize) %>% rename(Size="Freq") %>% mutate(EarlyP=EarlyP/Total,Lym=Lym/Total,ME=ME/Total,MK=MK/Total,Mye=Mye/Total,Total.norm=Total/Size,Total.norm_NPadj=Total/Stem) 
datatoplot<-datatoplot[order(datatoplot$Total.norm,decreasing=T),] %>% mutate(npClone=factor(npClone,levels=npClone))   
HSC_redeemR@CellMeta<-data.frame(VN=rowSums(HSC_redeemR@Cts.Mtx.bi)) %>% tibble::rownames_to_column("Cell")  %>% merge(HSC_redeemR@CellMeta,.)
depth_V<-HSC_redeemR@CellMeta %>% group_by(Clone_merge) %>% dplyr::summarise(Depth=mean(meanCov),VN=mean(VN)) %>% mutate(npClone=as.character(Clone_merge))
datatoplot<-merge(datatoplot,depth_V) %>% mutate(StemOdds=Stem/Total)
datatoplot.scale<-datatoplot %>% mutate(Size=Size,Lym=as.numeric(scale(Lym)),Mye=as.numeric(scale(Mye)),ME=as.numeric(scale(ME)),MK=as.numeric(scale(MK)), OutLevel.scale=as.numeric(scale(Total.norm)), OutLevel_NPadj.scale=as.numeric(scale(Total.norm_NPadj))) %>% 
dplyr::select(npClone,Lym,Mye,ME,MK,OutLevel.scale,OutLevel_NPadj.scale) %>% rename(Clone_merge="npClone")
return(list(Accuracy_HSC.BenchMark=Accuracy_HSC.BenchMark,ALLmeta.npClone=ALLmeta.npClone,ALLmeta.npClone.filtered=ALLmeta.npClone.filtered,datatoplot=datatoplot,datatoplot.scale=datatoplot.scale,np_mat=np_mat))
}

#' MakeDF4Regress  
#' Define a function to make two dataframe for regression analysis
#' This function was developed based on HSC_multiome_Het_2.ipynb
#' @param multiome_wrapper This outject should includes all and more than HSCs cells in redeemR
#' @param redeemR  scredeemR object for HSC
#' @param progeny_np run via ProgenyMapping_np
#' @param assay SCT for expression, ATAC for ATAC
#' @param useNPimputation default is T, use all cells called by network propagation, inaddition to the top cells in redeemR
#' @param maxcloneUMI default is 10, Only include genes, in the max clone the expression greater than 10
#' @return list(mtx.clone=mtx.clone,mtx.clone.norm.scale=mtx.clone.norm.scale)
#' @export
MakeDF4Regress<-function(multiome_wrapper=Donor04_HSC_Multiome_wrapper,
                         redeemR=DN4_stemcell_redeemR.seed.sensitive,
                         progeny_np=DN4_HSC_LSI_progeny_np,
                         assay="SCT",useNPimputation=T,maxcloneUMI=10){
## Get mtx
mtx<-multiome_wrapper@assays[[assay]]@counts %>% t %>% as.matrix
row.names(mtx)<-gsub("Cell","",row.names(mtx))
## Compute CloneInfo (Cell|npClone)
if(useNPimputation){
    CloneInfo<-subset(progeny_np$ALLmeta.npClone,STD.CellType=="HSC") %>% .[,c("Cell","npClone")] %>% mutate(Cell_RNA=Translate_simple_ATAC2RNA(Cell)) 
}else{
    CloneInfo<-redeemR@CellMeta[,c("Cell","Clone_merge")]%>% mutate(Cell_RNA=Translate_simple_ATAC2RNA(Cell)) %>% rename(npClone="Clone_merge")   
}
CloneInfo<-CloneInfo %>%tibble::remove_rownames() %>% tibble::column_to_rownames("Cell_RNA") %>% select(-Cell)

## Compute mtx.clone.norm(npClone|gene1|gene2|...)
mtx.clone<-Tomerge_v2(mtx,CloneInfo) %>% group_by(npClone) %>% dplyr::summarise_all(sum)
mtx.clone<-mtx.clone[!is.na(mtx.clone$npClone),]  ## Becasue some cells is NA when assigning clones
clonesize.umi<-mtx.clone[,-1] %>% rowSums ## The order of the clonesize.umi is consistent with the rows of mtx.clone 
mtx.clone.norm<-1e6*mtx.clone[,-1]/clonesize.umi 
mtx.clone.norm.scale<-cbind(npClone=mtx.clone[,1],scale(mtx.clone.norm)) 

## Only include genes with at least maxcloneUMI in max clone
mtx.clone[,-1] %>% apply(.,2,max) -> maxcloneExpr
features<-names(maxcloneExpr)[maxcloneExpr>=maxcloneUMI]
mtx.clone<-mtx.clone[,c("npClone",features)] %>% mutate(npClone=as.character(npClone)) %>% cbind(TotalUMI=clonesize.umi,.)
mtx.clone.norm.scale<-mtx.clone.norm.scale[,c("npClone",features)] %>% mutate(npClone=as.character(npClone))

## Merge in the lineage bias and output information
Lin_Out<-progeny_np$datatoplot.scale %>% rename(npClone="Clone_merge")
mtx.clone<-merge(Lin_Out,mtx.clone)
mtx.clone.norm.scale<-merge(Lin_Out,mtx.clone.norm.scale)
return(list(mtx.clone=mtx.clone,mtx.clone.norm.scale=mtx.clone.norm.scale))
}


#' Run_Lin_regression
#' 
#' Firstly used in HSC_multiome_Het_2.ipynb
#' @param LinOut  produced by MakeDF4Regress
#' @param regress_factor default is c("OutLevel.scale","OutLevel_NPadj.scale","Lym","Mye","MK","ME")
#' @param n.cores  default is 8
#' @export
#' @import foreach doParallel
Run_Lin_regression<-function(LinOut,regress_factor=c("OutLevel.scale","OutLevel_NPadj.scale","Lym","Mye","MK","ME"),n.cores=8){
require(dplyr)
require(foreach)
require(doParallel)
require(qvalue)
colnames(LinOut$mtx.clone.norm.scale)<-gsub("-","_",colnames(LinOut$mtx.clone.norm.scale))
colnames(LinOut$mtx.clone.norm.scale)<-gsub("/","_",colnames(LinOut$mtx.clone.norm.scale))
genes=colnames(LinOut$mtx.clone.norm.scale)[8:ncol(LinOut$mtx.clone.norm.scale)]
my.cluster <- parallel::makeCluster(n.cores)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
res<-foreach(gene=genes[1:length(genes)]) %dopar%{
    ps<-c()
    slopes<-c()
    for(i in 1:length(regress_factor)){
        f<-as.formula(paste(gene,"~",regress_factor)[i])
        md<-lm(f,data=LinOut$mtx.clone.norm.scale)
        md.summary<-summary(md)
        slope<-md.summary$coefficients[2,1]
        p<-md.summary$coefficients[2,4]
        slopes<-c(slopes,slope)
        ps<-c(ps,p)
    }
    return(list(Gene=gene,ps=ps,slopes=slopes))
}
parallel::stopCluster(cl = my.cluster)
genes<-sapply(res,function(x){x[[1]]})
ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame
slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame
genes<-gsub("_","-",genes)
row.names(ps)<-genes
row.names(slopes)<-genes
colnames(ps)<-regress_factor
colnames(slopes)<-regress_factor
qs<-apply(ps,2,function(x){qvalue(x)$qvalues})  %>% as.data.frame  
return(list(ps=ps,qs=qs,slopes=slopes))
}

#' Run_Lin_regression_poi
#' Firstly used in HSC_multiome_Het_2.ipynb
#' This function was developed based on 
#' @param LinOut produced by MakeDF4Regress
#' @param regress_factor  default is c("OutLevel.scale","OutLevel_NPadj.scale","Lym","Mye","MK","ME")
#' @param n.cores  =8
#' @export
#' @import foreach doParallel 
Run_Lin_regression_poi<-function(LinOut,regress_factor=c("OutLevel.scale","OutLevel_NPadj.scale","Lym","Mye","MK","ME"),n.cores=8){
require(dplyr)
require(foreach)
require(doParallel)
require(qvalue)
colnames(LinOut$mtx.clone)<-gsub("-","_",colnames(LinOut$mtx.clone))
colnames(LinOut$mtx.clone)<-gsub("/","_",colnames(LinOut$mtx.clone))
genes=colnames(LinOut$mtx.clone)[9:ncol(LinOut$mtx.clone)]
my.cluster <- parallel::makeCluster(n.cores)
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
res<-foreach(gene=genes[1:length(genes)]) %dopar%{
    ps<-c()
    slopes<-c()
    for(i in 1:length(regress_factor)){
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(TotalUMI)")[i])
        md<-glm(f,data=LinOut$mtx.clone,family=poisson(link="log"))
        md.summary<-summary(md)
        slope<-md.summary$coefficients[2,1]
        p<-md.summary$coefficients[2,4]
        slopes<-c(slopes,slope)
        ps<-c(ps,p)
    }
    return(list(Gene=gene,ps=ps,slopes=slopes))
}
parallel::stopCluster(cl = my.cluster)
genes<-sapply(res,function(x){x[[1]]})
ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame
slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame
genes<-gsub("_","-",genes)
row.names(ps)<-genes
row.names(slopes)<-genes
colnames(ps)<-regress_factor
colnames(slopes)<-regress_factor
qs<-apply(ps,2,function(x){qvalue(x)$qvalues})  %>% as.data.frame   
return(list(ps=ps,qs=qs,slopes=slopes))
}


#########################################Three plotting functions for HSC-progny analysis


#' plot_npSummary to assess the outputlevel
#'
#' @param npresult from ProgenyMapping_np
#' @param orderby Normalize by, so far can work with "Total.norm" and "Total.norm_NPadj"
#' @param pre Any short description for this plot to print with the plot
#' @export
## Define plotting functions
plot_npSummary<-function (npresult,orderby="Total.norm",pre)  # or Total.norm_NPadj
{
npsummary=npresult$datatoplot
npsummary<-npsummary[order(npsummary[,orderby],decreasing=T),] %>% mutate(npClone=factor(npClone,levels=npClone))    
p1 <- ggplot(npsummary) + aes(Size, Total) + geom_point() + ggtitle("ClonalSize VS Progeny number") + theme_bw()
p2 <- ggplot(npsummary) + aes_string("Size", orderby) + geom_point() + ggtitle("ClonalSize VS Normaized Progeny number") + theme_bw()
p.Size_VS_Stem<-ggplot(npsummary)+aes(Size,Stem)+geom_point()+theme_bw()+ ggtitle("ClonalSize VS AssignClonalSize")
p3 <- ggplot(npsummary) + aes_string("npClone", orderby, fill = "Depth") + geom_bar(stat = "identity") + theme_bw()+scale_fill_gradient2(low="grey",mid="white",high="red")+theme(axis.title.x = element_blank())
p4 <- ggplot(npsummary) + aes(npClone, "", color = Size, size = Size) + geom_point() + geom_point() + theme_bw() + scale_color_gradient(low = "darkgreen", high = "red") + ggtitle("CloneSize")+theme(axis.title.x = element_blank())
p4.assign <- ggplot(npsummary) + aes(npClone, "", color = Stem, size = Stem) + geom_point() + geom_point() + theme_bw() + scale_color_gradient(low = "darkgreen", high = "red") + ggtitle("AssignCloneSize")+theme(axis.title.x = element_blank())
p5 <- ggplot(npsummary) + aes(npClone, "", color = VN, size = VN) + geom_point() + geom_point() + theme_bw() + scale_color_viridis() + ggtitle("Variant #")+theme(axis.title.x = element_blank())
#p6 <- ggplot(npsummary) + aes(npClone, "", color = StemOdds, size = StemOdds) + geom_point() + geom_point() + theme_bw() + scale_color_viridis() + ggtitle("Stem odds")
l<-rbind(c(1,2,3),c(1,2,3),c(1,2,3),rep(4,3),rep(4,3),rep(4,3),rep(4,3),rep(4,3),rep(4,3),rep(5,3),rep(6,3),rep(7,3))
if(orderby=="Total.norm"){
    t<-"Normalized by Clone size"
}else if(orderby=="Total.norm_NPadj"){
    t<-"Normalized by Assigned Clonal Size"
}
grid.arrange(p1,p2,p.Size_VS_Stem,p3,p4,p4.assign,p5,layout_matrix=l,top=paste(pre,t))  
}

#' plot_npSummary to assess the outputlevel vs lineage bias, normalize by assigned
#'
#' @param datatoplot   A slot from the result of ProgenyMapping_np : datatoplot.scale
#' @param pre Any short description for this plot to print with the plot
#' @export
## Define plotting functions
Runplot_scale_2<-function(datatoplot=DN4_HSC_LSI_progeny$LineageSummary$output_lineage.summary.pct.scale,pre){
require(ggrepel)
require(ggpubr)
    p1<-datatoplot %>% mutate(LB=ifelse(abs(OutLevel_NPadj.scale)>1.28 | abs(Mye)>1.28,as.character(Clone_merge),"")) %>%
ggplot()+aes(OutLevel_NPadj.scale,Mye,label=LB)+geom_point()+geom_vline(xintercept = c(-1.28,1.28),linetype=2)+geom_hline(yintercept = c(-1.28,1.28),linetype=2)+
geom_text_repel()+theme_classic()+geom_smooth(method="lm")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 1)#+stat_regline_equation(aes(label = ..p.label..))#+stat_regline_equation(aes(label = ..rr.label..))

p2<-datatoplot %>% mutate(LB=ifelse(abs(OutLevel_NPadj.scale)>1.28 | abs(Lym)>1.28,as.character(Clone_merge),"")) %>%
ggplot()+aes(OutLevel_NPadj.scale,Lym,label=LB)+geom_point()+geom_vline(xintercept = c(-1.28,1.28),linetype=2)+geom_hline(yintercept = c(-1.28,1.28),linetype=2)+
geom_text_repel()+theme_classic()+geom_smooth(method="lm")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 1)#+stat_regline_equation(aes(label = ..rr.label..))

p3<-datatoplot %>% mutate(LB=ifelse(abs(OutLevel_NPadj.scale)>1.28 | abs(MK)>1.28,as.character(Clone_merge),"")) %>%
ggplot()+aes(OutLevel_NPadj.scale,MK,label=LB)+geom_point()+geom_vline(xintercept = c(-1.28,1.28),linetype=2)+geom_hline(yintercept = c(-1.28,1.28),linetype=2)+
geom_text_repel()+theme_classic()+geom_smooth(method="lm")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 1)#+stat_regline_equation(aes(label = ..rr.label..))

p4<-datatoplot %>% mutate(LB=ifelse(abs(OutLevel_NPadj.scale)>1.28 | abs(ME)>1.28,as.character(Clone_merge),"")) %>%
ggplot()+aes(OutLevel_NPadj.scale,ME,label=LB)+geom_point()+geom_vline(xintercept = c(-1.28,1.28),linetype=2)+geom_hline(yintercept = c(-1.28,1.28),linetype=2)+
geom_text_repel()+theme_classic()+geom_smooth(method="lm")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 1)#+stat_regline_equation(aes(label = ..r.label..))#+stat_regline_equation(aes(label = ..rr.label..))
grid.arrange(p1,p2,p3,p4,top=paste(pre,"Normalize by assigned size"))
}

#' plot_npSummary to assess the outputlevel vs lineage bias, normalize by HSC original clone size
#'
#' @param datatoplot   A slot from the result of ProgenyMapping_np : datatoplot.scale
#' @param pre Any short description for this plot to print with the plot
#' @export
## Define plotting functions
Runplot_scale_3<-function(datatoplot=DN4_HSC_LSI_progeny$LineageSummary$output_lineage.summary.pct.scale,pre){
require(ggrepel)
require(ggpubr)
    p1<-datatoplot %>% mutate(LB=ifelse(abs(OutLevel.scale)>1.28 | abs(Mye)>1.28,as.character(Clone_merge),"")) %>%
ggplot()+aes(OutLevel.scale,Mye,label=LB)+geom_point()+geom_vline(xintercept = c(-1.28,1.28),linetype=2)+geom_hline(yintercept = c(-1.28,1.28),linetype=2)+
geom_text_repel()+theme_classic()+geom_smooth(method="lm")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 1)#+stat_regline_equation(aes(label = ..p.label..))#+stat_regline_equation(aes(label = ..rr.label..))

p2<-datatoplot %>% mutate(LB=ifelse(abs(OutLevel.scale)>1.28 | abs(Lym)>1.28,as.character(Clone_merge),"")) %>%
ggplot()+aes(OutLevel.scale,Lym,label=LB)+geom_point()+geom_vline(xintercept = c(-1.28,1.28),linetype=2)+geom_hline(yintercept = c(-1.28,1.28),linetype=2)+
geom_text_repel()+theme_classic()+geom_smooth(method="lm")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 1)#+stat_regline_equation(aes(label = ..rr.label..))

p3<-datatoplot %>% mutate(LB=ifelse(abs(OutLevel.scale)>1.28 | abs(MK)>1.28,as.character(Clone_merge),"")) %>%
ggplot()+aes(OutLevel.scale,MK,label=LB)+geom_point()+geom_vline(xintercept = c(-1.28,1.28),linetype=2)+geom_hline(yintercept = c(-1.28,1.28),linetype=2)+
geom_text_repel()+theme_classic()+geom_smooth(method="lm")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 1)#+stat_regline_equation(aes(label = ..rr.label..))

p4<-datatoplot %>% mutate(LB=ifelse(abs(OutLevel.scale)>1.28 | abs(ME)>1.28,as.character(Clone_merge),"")) %>%
ggplot()+aes(OutLevel.scale,ME,label=LB)+geom_point()+geom_vline(xintercept = c(-1.28,1.28),linetype=2)+geom_hline(yintercept = c(-1.28,1.28),linetype=2)+
geom_text_repel()+theme_classic()+geom_smooth(method="lm")+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),label.x = 1)#+stat_regline_equation(aes(label = ..r.label..))#+stat_regline_equation(aes(label = ..rr.label..))
grid.arrange(p1,p2,p3,p4,top=paste(pre,"Normalized by Size"))
}


#' plot_npSummary to plot the lineage composition
#'
#' @param npresult from ProgenyMapping_np
#' @param pre Any short description for this plot to print with the plot
#' @export
## Define plotting functions
LineageBiasPlot<-function(npresult,pre){
library("wesanderson")
library(RColorBrewer)
library(ggsci)
datatoplot=npresult$datatoplot
datatoplot<-datatoplot %>% dplyr::select(npClone,MK,Lym,ME,Total.norm_NPadj,Mye) %>%
mutate(Lym=Lym/(Lym+Mye+ME+MK),Mye=Mye/(Lym+Mye+ME+MK),ME=ME/(Lym+Mye+ME+MK),MK=MK/(Lym+Mye+ME+MK))
ps<-list()
for(orderby in c("MK","Lym","ME","Mye")){
p<-datatoplot[order(datatoplot[,orderby],decreasing=F),] %>% mutate(npClone=factor(npClone,levels=as.character(npClone))) %>% dplyr::select(-Total.norm_NPadj)%>% reshape2::melt() %>%
ggplot()+aes(npClone,value,fill=variable)+geom_bar(stat="identity",position="fill",color="black",size=0.2)+scale_fill_npg()+ggtitle(paste("Order by:",orderby))+theme_pubr()+theme(axis.text.x  = element_text(angle = 45))
ps<-c(ps,list(p)) 
}
options(repr.plot.width=20, repr.plot.height=12,repr.plot.res=100)
grid.arrange(grobs=ps,top=pre)
    
}
