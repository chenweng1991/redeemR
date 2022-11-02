#' Wrap Seurat RNA clustering
#'
#' This function allows you to perform standard sc-RNA clustering
#' @param mtx sparse Matrix of class "dgCMatrix", each row is a gene, each column is a cell,
#' @param exp The name of this sample/experiment
#' @param res clustering resolution, default=0.5
#' @return this returns seurat object with RNA clustering
#' @examples
#' bmmc.data=Read10X(data.dir = "XX/CellRanger/Donor01_BMMC_1/outs/filtered_feature_bc_matrix")
#' docluster_GEM(mtx=bmmc.data$`Gene Expression`,exp="DN1_BMMC1")
#' @export
GEM_Wrapper<-function(mtx=bmmc.data$`Gene Expression`,exp="DN1_BMMC1",res=0.5){
require(Seurat)
ob <- CreateSeuratObject(counts = mtx, project = exp, min.cells = 3, min.features = 200)
ob <- NormalizeData(ob, normalization.method = "LogNormalize", scale.factor = 10000)
ob <- FindVariableFeatures(ob, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ob)
ob <- ScaleData(ob, features = all.genes)
ob <- RunPCA(ob, features = VariableFeatures(object = ob))
ob <- FindNeighbors(ob, dims = 1:10)
ob <- FindClusters(ob, resolution = res)
ob <- RunUMAP(ob, dims = 1:10)
return(ob)
}

#' Wrap Seurat ATAC clustering
#'
#' This function allows you to perform standard sc-ATAC clustering
#' @param MTX sparse Matrix of class "dgCMatrix", each row is a peak, each column is a cell,
#' @param res clustering resolution, default=0.5
#' @return this returns seurat object with ATAC clustering
#' @examples
#' bmmc.filtered.atac<-SeuratLSIClustering(PeakVSCell.filtered.Mtx) #each row is a peak, each column is a cell,
#' @export
ATAC_Wrapper<-function(MTX,res=0.3,dim1=1, dim2=20){
require(Signac)
Cell_Variant.seurat<-CreateSeuratObject(counts = MTX, assay = "ATAC")
VariableFeatures(Cell_Variant.seurat) <- row.names(Cell_Variant.seurat) #names(which(Matrix::rowSums(Cell_Variant.seurat) > 100))
Cell_Variant.seurat <- RunTFIDF(Cell_Variant.seurat, n = 50)
Cell_Variant.seurat<- FindTopFeatures(Cell_Variant.seurat, min.cutoff = 'q0')
Cell_Variant.seurat <- RunSVD(Cell_Variant.seurat, n = 50)
Cell_Variant.seurat <- RunUMAP(Cell_Variant.seurat, reduction = "lsi", dims = dim1:dim2)
Cell_Variant.seurat <- FindNeighbors(Cell_Variant.seurat,reduction ="lsi"  ,dims = dim1:dim2)
Cell_Variant.seurat <- FindClusters(Cell_Variant.seurat, resolution = res)
}



#' Wrap Seurat Multiomics clustering
#'
#' This function allows you to perform standard sc-multiome clustering
#' @param path this should be the path to the cell-ranger results XX/outs
#' @param atacmin minimum atac fragment for each cell, default is 1000
#' @param umimin minimum rna umi for each cell, default is 1000
#' @param cellID to be used for input(useful for re-clustering), default is NULL which will use the info from path/per_barcode_metrics.csv
#' @return this returns seurat object with both RNA and ATAC
#' @examples
#' Multi_Wrapper(path="XX/CellRanger/Donor01_BMMC_1/outs/")
#' @export
Multi_Wrapper<-function(path="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/CellRanger/Donor01_CD34_1/outs",atacmin=1000,umimin=1000,CellID=NULL){
require(Seurat)
require(Signac)
require(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
require(dplyr)
require(ggplot2)
inputdata.10x <- Read10X_h5(paste(path,"/raw_feature_bc_matrix.h5",sep=""))
# per_barcode_metrics<-read.csv(paste(path,"/per_barcode_metrics.csv",sep=""))
if(length(CellID)==0){
CellID<-read.table(paste(path,"/filtered_feature_bc_matrix/barcodes.tsv.gz",sep=""))$V1
}
# Extract rna and atac counts
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Only use standard chromasome for atac counts
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- as.character(seqnames(grange.counts)) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
# Filter atac and rna counts
rna_counts.filtered<-rna_counts[,CellID]
atac_counts.filtered<-atac_counts[,CellID]
# Use RNA to create the default object
ob<-CreateSeuratObject(counts = rna_counts.filtered)
# Create chrome_assay
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations)<-"UCSC"
genome(annotations) <- "hg38"
frag.file <- paste(path,"/atac_fragments.tsv.gz",sep="")
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts.filtered,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
)
# Add chrom_assay
ob[["ATAC"]]<-chrom_assay
## Further filter the object
ob <- subset(
  x = ob,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > atacmin &
    nCount_RNA < 25000 &
    nCount_RNA > umimin
)
# RNA analysis
DefaultAssay(ob) <- "RNA"
ob <- SCTransform(ob, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(ob) <- "ATAC"
ob <- RunTFIDF(ob)
ob <- FindTopFeatures(ob, min.cutoff = 'q0')
ob <- RunSVD(ob)
ob <- RunUMAP(ob, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
ob <- FindMultiModalNeighbors(ob, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
ob <- RunUMAP(ob, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ob <- FindClusters(ob, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
return(ob)
}


