#' Function to summarize the depth (Total that passed Q30)
#'
#' This function allows you to summarize the depth
#' @param path The XX/final folder, the output from mitoV
#' @param Processed Boolean variable(Default T), if true directly readRDS("depth.RDS") or, generate and saveout "depth.RDS"
#' @return this returns depth which is a list of 4 list(Total/VerySensitive/Sensitive/Specific), each contains 2 df, summarize mito coverage by Pos/Cell
#' @examples
#' WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/MTenrichCombine/Enrich/CW_mgatk/final"
#' DN1CD34_1.depth<-DepthSummary(WD,Processed = T)
#' @export
#' @import dplyr
DepthSummary<-function(path,Processed=T){
setwd(path)
if(Processed){
    depth<-readRDS("depth.RDS")
}else{
    QualifiedTotalCts<-read.table("QualifiedTotalCts")
    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V3))
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V3))
    depth_Total<-list(Pos.MeanCov,Cell.MeanCov)

    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V4))
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V4))
    depth_VerySensitive<-list(Pos.MeanCov,Cell.MeanCov)

    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V5))
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V5))
    depth_Sensitive<-list(Pos.MeanCov,Cell.MeanCov)

    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V6))
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V6))
    depth_Specific<-list(Pos.MeanCov,Cell.MeanCov)

    depth<-list(Total=depth_Total,VerySensitive=depth_VerySensitive,Sensitive=depth_Sensitive,Specific=depth_Specific)
    saveRDS(depth,"depth.RDS")
return(depth)
}
}


#' Function to generate GTS summary
#'
#' This function allows you to summarize the meta data for each genotyped variant
#' @param RawGenotypes Well-named "RawGenotypes.Sensitive.StrandBalance" file in function CW_mgatk.read
#' @param filterN Boolean variable, if true filter out the variant with "N"
#' @return Genotypes.summary a dataframe that summarize several metrics for each genotype
#' @examples
#' Usually used inside of function CW_mgatk.read
#' @export
#' @import dplyr
GTSummary<-function(RawGenotypes,filterN=T){ ## At this moment, the context with N is probably prone to error due to mapping, in the future should work on realignment
# Make Depth dictionary
Depth<-unique(RawGenotypes[,c("Cell","Pos","Depth")])
Depthdic<-Depth$Depth
names(Depthdic)<-paste(Depth$Cell, Depth$Pos,sep="")
# Summarise
Genotypes.summary<-table(paste(RawGenotypes$Cell,RawGenotypes$Variants,sep="_")) %>% as.data.frame()
Genotypes.summary$Cell<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){x[1]})
Genotypes.summary$Variants<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){paste(x[2:4],collapse="_")})
Genotypes.summary$cellPos<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){paste(x[1:2],collapse="")})
Genotypes.summary$depth<-Depthdic[Genotypes.summary$cellPos]
Genotypes.summary<-Genotypes.summary[,c("Var1","Cell","Variants","Freq","depth")]
Genotypes.summary$Type<-strsplit(Genotypes.summary$Variants,"_") %>% sapply(.,function(x){paste(x[2],x[3],sep="_")})
Genotypes.summary$Context<-ContextsDic[strsplit(Genotypes.summary$Variants,"_") %>% sapply(.,function(x){x[1]})]
if(filterN){
    Genotypes.summary<-subset(Genotypes.summary,!grepl("N",Genotypes.summary$Context))
}
return(Genotypes.summary)
}




#' Function to read in mitoV outputs
#'
#' This function allows you to read raw data from XX/final folder, the output from mitoV
#' @param path The XX/final folder, the output from mitoV
#' @param Processed Boolean variable (Default F), if true directly readRDS("VariantsGTSummary.RDS") or, generate and saveout "VariantsGTSummary.RDS"
#' @return this returns depth which is a list of 4 df (Total/VerySensitive/Sensitive/Specific), each is a genotype summary
#' @examples
#' WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/MTenrichCombine/Enrich/CW_mgatk/final"
#' DN1CD34_1.VariantsGTSummary<-CW_mgatk.read(WD,Processed =T)
#' @export
CW_mgatk.read<-function(path,Processed=F){
setwd(path)
if(Processed){
    VariantsGTSummary<-readRDS("VariantsGTSummary.RDS")
}else{
    RawGenotypes.Sensitive.StrandBalance<-read.table("RawGenotypes.Sensitive.StrandBalance")
    RawGenotypes.VerySensitive.StrandBalance<-read.table("RawGenotypes.VerySensitive.StrandBalance")
    RawGenotypes.Specific.StrandBalance<-read.table("RawGenotypes.Specific.StrandBalance")
    RawGenotypes.Total.StrandBalance<-read.table("RawGenotypes.Total.StrandBalance")
    GiveName<-c("UMI","Cell","Pos","Variants","Call","Ref","FamSize","GT_Cts","CSS","DB_Cts","SG_Cts","Plus","Minus","Depth")
    colnames(RawGenotypes.Sensitive.StrandBalance)<-GiveName
    colnames(RawGenotypes.VerySensitive.StrandBalance)<-GiveName
    colnames(RawGenotypes.Specific.StrandBalance)<-GiveName
    colnames(RawGenotypes.Total.StrandBalance)<-GiveName
    Specific.GTSummary<-GTSummary(RawGenotypes.Specific.StrandBalance)
    Sensitive.GTSummary<-GTSummary(RawGenotypes.Sensitive.StrandBalance)
    VerySensitive.GTSummary<-GTSummary(RawGenotypes.VerySensitive.StrandBalance)
    Total.GTSummary<-GTSummary(RawGenotypes.Total.StrandBalance)
    ##Calculate heteroplasmy
    Specific.GTSummary$hetero<-with(Specific.GTSummary,Freq/depth)
    Sensitive.GTSummary$hetero<-with(Sensitive.GTSummary,Freq/depth)
    VerySensitive.GTSummary$hetero<-with(VerySensitive.GTSummary,Freq/depth)
    Total.GTSummary$hetero<-with(Total.GTSummary,Freq/depth)
    VariantsGTSummary<-list(Total=Total.GTSummary,VerySensitive=VerySensitive.GTSummary,Sensitive=Sensitive.GTSummary,Specific=Specific.GTSummary)
    saveRDS(VariantsGTSummary,"VariantsGTSummary.RDS")
}
return(VariantsGTSummary)
}

#' Function to filter variants
#'
#' This function allows you to filter variants
#' @param InputSummary The GTSummary file read in by function CW_mgatk.read
#' @param depth The .depth file by function DepthSummary
#' @param Rmvhomo Boolean (Default F) If true, remove the homozygous variants
#' @param Min_Cells Default 2, A qualified variant needs the minimum number of cells that have this variant
#' @param Max_Count_perCell Default 2,  A qualified variant needs to show at least 2 counts in one cell
#' @param QualifyCellCut Default 10, Minimum depth for a qualified cell
#' @return this returns feature.list
#' @examples
#' plot_variant(DN1CD34_1.VariantsGTSummary,DN1CD34_1.Variants.feature.lst,depth=DN1CD34_1.depth,cat=
#' c("Total","VerySensitive","Sensitive","Specific"),p4xlim = 30)
#' @export
#' @import dplyr
Vfilter_v3<-function(InputSummary,depth,Rmvhomo=F,Min_Cells=2, Max_Count_perCell=2,QualifyCellCut=10){
    CV<-function(x){
    var(x)/mean(x)
    }
Names<-names(InputSummary)
feature.list<-list()
for(i in Names){
VariantFeature0<- InputSummary[[i]] %>% group_by(Variants) %>% dplyr::summarise(CellN=n(),PositiveMean=mean(hetero),maxcts=max(Freq),CV=CV(hetero),TotalVcount=sum(Freq))
VariantFeature0$pos<-strsplit(VariantFeature0$Variants,"_") %>% sapply(.,function(x){x[1]}) %>% as.numeric
VariantFeature0<-merge(VariantFeature0,depth[[i]][[1]],by.x="pos",by.y="V2")   ## This generate different meanCov for each threahold
VariantFeature0$TotalCov<-length(unique(InputSummary[[i]]$Cell))*VariantFeature0$meanCov
VariantFeature0$VAF<-VariantFeature0$TotalVcount/VariantFeature0$TotalCov
qualifiedCell<-subset(depth[["Total"]][[2]],meanCov>=QualifyCellCut)[,1,drop=T]  ## Filter Qualified cell based on total depth
InputSummary.qualified<-subset(InputSummary[[i]],Cell %in% qualifiedCell)
VariantFeature<- InputSummary.qualified %>% group_by(Variants) %>% dplyr::summarise(CellN=n(),PositiveMean=mean(hetero),maxcts=max(Freq),CV=CV(hetero),TotalVcount=sum(Freq))
print(paste(i,":\n",nrow(VariantFeature0),"variants to start"))
print(paste(nrow(VariantFeature),"variants after remove low quality cells"))
VariantFeature$CellNPCT<-VariantFeature$CellN/length(unique(InputSummary.qualified$Cell))
VariantFeature<-merge(VariantFeature[,c("Variants","CellN","PositiveMean","maxcts","CellNPCT")],VariantFeature0[,c("Variants","TotalVcount","TotalCov","VAF","CV")],by="Variants")
#HomoVariants<-subset(VariantFeature,VAF>0.9 & CV<0.01)$Variants
HomoVariants<-subset(VariantFeature,CellNPCT>0.75 & PositiveMean>0.75 & CV<0.01)$Variants
VariantFeature$HomoTag<-ifelse(VariantFeature$Variants %in% HomoVariants,"Homo","Hetero")
if (Rmvhomo){
    VariantFeature<-subset(VariantFeature,!Variants %in% HomoVariants)
    print(paste(length(HomoVariants),"Homoplasmy variants to remove"))
    print(HomoVariants)
}else{
    print(paste("Tag Homoplasmy:",HomoVariants))
}
out<-subset(VariantFeature,CellN>=Min_Cells & maxcts>=Max_Count_perCell)
print(paste("After filtering,",nrow(out), "Variants left"))
print("\n\n")
feature.list<-c(feature.list,list(out))
}
    names(feature.list)<-Names
    return(feature.list)
}
