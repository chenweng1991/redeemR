#' Function to summarize the depth (Total that passed Q30)
#'
#' This function allows you to summarize the depth
#' @param path The XX/final folder, the output from mitoV
#' @param Processed Boolean variable(Default T), if true directly readRDS("depth.RDS") or, generate and saveout "depth.RDS"
#' @param CellSubset A vector of ATAC cell names for subsetting, default is NA
#' @param only_Total Default is T, Only return total depth summary. Don't care about depth in different quality filtering
#' @return this returns depth which is a list of 4 list(Total/VerySensitive/Sensitive/Specific), each contains 2 df, summarize mito coverage by Pos/Cell
#' @examples
#' WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/MTenrichCombine/Enrich/CW_mgatk/final"
#' DN1CD34_1.depth<-DepthSummary(WD,Processed = T)
#' @export
#' @import dplyr
DepthSummary<-function(path,CellSubset=NA,only_Total=T){
    message("By default only total depth is summarized")
    message("deprecated")
    QualifiedTotalCts<-read.table(paste(path,"/QualifiedTotalCts",sep=""))
    if(!is.na(CellSubset)){
      print("Will subset cells...")
      print(paste(length(unique(QualifiedTotalCts$V1)),"Cells in QualifiedTotalCts"))
      print(paste(length(CellSubset),"Cells in the provided subset"))
      print(paste(length(which(CellSubset %in% unique(QualifiedTotalCts$V1))),"Subset cells overlap in QualifiedTotalCts"))
      QualifiedTotalCts<-subset(QualifiedTotalCts,V1 %in% CellSubset)
    }else{
      print("Use all cells")
      print(paste(length(unique(QualifiedTotalCts$V1)),"Cells in QualifiedTotalCts"))
    }
    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V3))
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V3))
    depth_Total<-list(Pos.MeanCov=Pos.MeanCov,Cell.MeanCov=Cell.MeanCov)
    if (only_Total){
        return(depth_Total)
    }else{
    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V4))
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V4))
    depth_LS<-list(Pos.MeanCov=Pos.MeanCov,Cell.MeanCov=Cell.MeanCov)

    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V5))
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V5))
    depth_S<-list(Pos.MeanCov=Pos.MeanCov,Cell.MeanCov=Cell.MeanCov)

    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V6))
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V6))
    depth_VS<-list(Pos.MeanCov=Pos.MeanCov,Cell.MeanCov=Cell.MeanCov)

    depth<-list(Total=depth_Total,VerySensitive=depth_LS,Sensitive=depth_S,Specific=depth_VS)
    return(depth)
    }
}



#' Function to generate GTS summary
#'
#' This function allows you to summarize the meta data for each genotyped variant
#' @param RawGenotypes Well-named "RawGenotypes.Sensitive.StrandBalance" file in function redeemR.read or CW_mgatk.read
#' @param filterN Boolean variable, if true filter out the variant with "N"
#' @return Genotypes.summary a dataframe that summarize several metrics for each genotype
#' @examples Usually used inside of function CW_mgatk.read
#' @export
#' @import dplyr
GTSummary<-function(RawGenotypes,filterN=T){ ## At this moment, the context with N is probably prone to error due to mapping, in the future should work on realignment
# Make Depth dictionary
data(ContextsDic)
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




#' Function to read in redeemV outputs
#'
#' This function allows you to read raw data from XX/final folder, the output from redeemV
#' It process the data same way as CW_mgatk.read but need to specify one threadhold(thr)
#' @param path The XX/final folder, the output from mitoV
#' @param thr The thredhold of filtering T(Total),LS(Less Stringent:c=0.75,a1=2,a2=1), S(Stringent:c=0.75,a1=3,a2=2), VS(Very Stringent:c=0.75,a1=4,a2=3)"
#' @param Processed Boolean variable (Default F), if true directly readRDS("VariantsGTSummary.RDS") or, generate and saveout "VariantsGTSummary.RDS"
#' @return this returns depth which is a list of 4 df (Total/VerySensitive/Sensitive/Specific), each is a genotype summary
#' @examples WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/MTenrichCombine/Enrich/CW_mgatk/final"
#' DN1CD34_1.VariantsGTSummary<-CW_mgatk.read(WD,Processed =T)
#' @export
redeemR.read<-function(path,thr="S",Processed=F,rdsname="/VariantsGTSummary.RDS"){
if(Processed){
    VariantsGTSummary<-readRDS(paste(path,"/VariantsGTSummary.RDS",sep=""))
}else{
    if(missing(path)|missing(thr)){
        message("redeemR.read(path,thr)")
        message("missing variable path or thr")
        message("path is a string to the redeemV result folder that contains RawGenotypes.XX")
        message("thr is one of T,LS,S,VS:")
        message("T(Total),LS(Less Stringent:c=0.75,a1=2,a2=1), S(Stringent:c=0.75,a1=3,a2=2), VS(Very Stringent:c=0.75,a1=4,a2=3)")
        message("Term from redeemV is deprecated: VerySensitive equal to Less Stringent, Sensitive equal to Stringent, Specific equal to Very Stringent")
        return(NULL)
    }
    GiveName<-c("UMI","Cell","Pos","Variants","Call","Ref","FamSize","GT_Cts","CSS","DB_Cts","SG_Cts","Plus","Minus","Depth")
    if(thr=="T"){
        RawGenotypes<-read.table(paste(path,"/RawGenotypes.Total.StrandBalance",sep=""))
    }else if(thr=="LS"){
        RawGenotypes<-read.table(paste(path,"/RawGenotypes.VerySensitive.StrandBalance",sep=""))
    }else if(thr=="S"){
        RawGenotypes<-read.table(paste(path,"/RawGenotypes.Sensitive.StrandBalance",sep=""))
    }else if(thr=="VS"){
        RawGenotypes<-read.table(paste(path,"/RawGenotypes.Specific.StrandBalance",sep=""))
    } 
    colnames(RawGenotypes)<-GiveName
    VariantsGTSummary<-GTSummary(RawGenotypes)
    ##Calculate heteroplasmy
    VariantsGTSummary$hetero<-with(VariantsGTSummary,Freq/depth)
    attr(VariantsGTSummary,"thr")<-thr
    attr(VariantsGTSummary,"depth")<-DepthSummary(path)
    attr(VariantsGTSummary,"path")<-path
    saveRDS(VariantsGTSummary,paste(path,rdsname,sep=""))
    return(VariantsGTSummary)
}
}

#' Internal CV 
#'
#' This function allows you to read raw data from XX/final folder, the output from mitoV
#' @param x input a vector of numeric values
CV<-function(x){
    var(x)/mean(x)
}


#' Old Function to read in redeemV outputs 
#' It process the data same way as redeemR.read but simultanously reading in all threadhold as a list
#' This function allows you to read raw data from XX/final folder, the output from redeemV
#' @param path The XX/final folder, the output from mitoV
#' @param Processed Boolean variable (Default F), if true directly readRDS("VariantsGTSummary.RDS") or, generate and saveout "VariantsGTSummary.RDS"
#' @return this returns depth which is a list of 4 df (Total/VerySensitive/Sensitive/Specific), each is a genotype summary
#' @examples WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/MTenrichCombine/Enrich/CW_mgatk/final"
#' DN1CD34_1.VariantsGTSummary<-CW_mgatk.read(WD,Processed =T)
#' @export
CW_mgatk.read<-function(path,Processed=F){
message("CW_mgatk.read is the old version (simultanously reading in all threadhold as a list), it has been deprecated, please use redeemR.read")
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

#' Function to filter variants, deprecated
#'
#' This function allows you to filter variants,deprecated, use Vfilter_v4 instead
#' @param InputSummary The GTSummary file read in by function CW_mgatk.read
#' @param depth The .depth file by function DepthSummary
#' @param Rmvhomo Boolean (Default F) If true, remove the homozygous variants
#' @param Min_Cells Default 2, A qualified variant needs the minimum number of cells that have this variant
#' @param Max_Count_perCell Default 2,  A qualified variant needs to show at least 2 counts in one cell
#' @param QualifyCellCut Default 10, Minimum depth for a qualified cell
#' @return this returns feature.list
#' @examples
#' DN1CD34_1.Variants.feature.lst<-Vfilter_v3(InputSummary=DN1CD34_1.VariantsGTSummary,depth=DN1CD34_1.depth)
#' @export
#' @import dplyr
Vfilter_v3<-function(InputSummary,depth,Rmvhomo=F,Min_Cells=2, Max_Count_perCell=2,QualifyCellCut=10){
message("Vfilter_v3 is the old version that works with CW_mgatk.read , it has been deprecated, please use Vfilter_v4 with redeemR.read")    
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

#' Function to filter variants, v4
#'
#' This function allows you to filter variants
#' @param InputSummary The GTSummary file read in by function redeemR.read
#' @param depth The .depth file by function DepthSummary
#' @param Rmvhomo Boolean (Default F) If true, remove the homozygous variants
#' @param Min_Cells Default 2, A qualified variant needs the minimum number of cells that have this variant
#' @param Max_Count_perCell Default 2,  A qualified variant needs to show at least 2 counts in one cell
#' @param QualifyCellCut Default 10, Minimum depth for a qualified cell
#' @return this returns feature.list
#' @export
#' @import dplyr
Vfilter_v4<-function(InputSummary=VariantsGTSummary,Min_Cells=2, Max_Count_perCell=2,QualifyCellCut=10){   
VariantFeature0<- InputSummary %>% group_by(Variants) %>% dplyr::summarise(CellN=n(),PositiveMean=mean(hetero),maxcts=max(Freq),CV=CV(hetero),TotalVcount=sum(Freq))
VariantFeature0$pos<-strsplit(VariantFeature0$Variants,"_") %>% sapply(.,function(x){x[1]}) %>% as.numeric
VariantFeature0<-merge(VariantFeature0,attr(InputSummary,"depth")[["Pos.MeanCov"]],by.x="pos",by.y="V2")   ## This generate different meanCov for each threahold
VariantFeature0$TotalCov<-length(unique(InputSummary$Cell))*VariantFeature0$meanCov
VariantFeature0$totalVAF<-VariantFeature0$TotalVcount/VariantFeature0$TotalCov
qualifiedCell<-subset(attr(InputSummary,"depth")[["Cell.MeanCov"]],meanCov>=QualifyCellCut)[,1,drop=T]  ## Filter Qualified cell based on total depth
InputSummary.qualified<-subset(InputSummary,Cell %in% qualifiedCell)
VariantFeature<- InputSummary.qualified %>% group_by(Variants) %>% dplyr::summarise(CellN=n(),PositiveMean=mean(hetero),maxcts=max(Freq),CV=CV(hetero),TotalVcount=sum(Freq))
print(paste(nrow(VariantFeature0),"variants to start"))
print(paste(nrow(VariantFeature),"variants after remove low quality cells"))
VariantFeature$CellNPCT<-VariantFeature$CellN/length(unique(InputSummary.qualified$Cell))
VariantFeature<-merge(VariantFeature[,c("Variants","CellN","PositiveMean","maxcts","CellNPCT")],VariantFeature0[,c("Variants","TotalVcount","TotalCov","totalVAF","CV")],by="Variants")
HomoVariants<-subset(VariantFeature,CellNPCT>0.75 & PositiveMean>0.75 & CV<0.01)$Variants
VariantFeature$HomoTag<-ifelse(VariantFeature$Variants %in% HomoVariants,"Homo","Hetero")
print(paste("Tag Homoplasmy:",HomoVariants))
VariantFeature.filtered<-subset(VariantFeature,CellN>=Min_Cells & maxcts>=Max_Count_perCell)
print(paste("After filtering,",nrow(VariantFeature.filtered), "Variants left"))
attr(VariantFeature.filtered,"HomoVariants")<-HomoVariants
attr(VariantFeature.filtered,"Filter.Cell")<-c(AllCellN=nrow(attr(InputSummary,"depth")[["Cell.MeanCov"]]),QualifiedCellN=length(qualifiedCell))
attr(VariantFeature.filtered,"Filter.V")<-c(VN_total=nrow(VariantFeature0),VN_rmvLowQualityCell=nrow(VariantFeature),VN_filter=nrow(VariantFeature.filtered),VN_filter_rmvHomo=nrow(subset(VariantFeature.filtered,HomoTag!="Homo")))
return(VariantFeature.filtered)
}


#' Function to compute the reject rate(The filtering rate in concensus variant calling)
#'
#' This function allows you to computae the filtering rate for each single cell
#' @param ob The redeemR object
#' @return  a modified ob with RejectRate added to @CellMeta
#' @export
ComputeRejectRate<-function(ob){
WD=ob@attr$path
Total<-read.table(paste(WD,"/RawGenotypes.Total.StrandBalance",sep=""))
if(ob@para["Threhold"]=="T"){
    RawGenotypes<-read.table(paste(ob@attr$path,"/RawGenotypes.Total.StrandBalance",sep=""))
}else if(ob@para["Threhold"]=="LS"){
    RawGenotypes<-read.table(paste(ob@attr$path,"/RawGenotypes.VerySensitive.StrandBalance",sep=""))
}else if(ob@para["Threhold"]=="S"){
    RawGenotypes<-read.table(paste(ob@attr$path,"/RawGenotypes.Sensitive.StrandBalance",sep=""))
}else if(ob@para["Threhold"]=="VS"){
    RawGenotypes<-read.table(paste(ob@attr$path,"/RawGenotypes.Specific.StrandBalance",sep=""))
}    
res<-Total %>% group_by(V2) %>% dplyr::summarise(Total=n())
res<-RawGenotypes %>% group_by(V2) %>% dplyr::summarise(RejectRate=n()) %>% merge(res,.,by="V2")
res<-res %>% mutate(RejectRate=RejectRate/Total) %>% select(-Total)
colnames(res)[1]<-"Cell"
ob@CellMeta<-merge(ob@CellMeta,res,by="Cell")
message("RejectRate has been added to @CellMeta")
return(ob)
}

#' Function to compute the reject rate(The filtering rate in concensus variant calling)
#'
#' This function allows you to computae the filtering rate for each single cell
#' @param WD The path to the work space usually  XXX/mitoV/final
#' @return a dataframe that store the percentage of variant in a given threahold again total
#' @examples
#' DN9_BMMC_RejectRate<-ComputeRejectRate("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9/Donor9/DN9_BMMC/MTenrichCombine/mitoV/final/")
#' @export
ComputeRejectRate_legacy<-function(WD){
message("This has been deprecated, please use ComputeRejectRate that takes redeemR object as input")
Total<-read.table(paste(WD,"/RawGenotypes.Total.StrandBalance",sep=""))
VerySensitive<-read.table(paste(WD,"/RawGenotypes.VerySensitive.StrandBalance",sep=""))
Sensitive<-read.table(paste(WD,"/RawGenotypes.Sensitive.StrandBalance",sep=""))
Specific<-read.table(paste(WD,"/RawGenotypes.Specific.StrandBalance",sep=""))
res<-Total %>% group_by(V2) %>% dplyr::summarise(Total=n())
res<-VerySensitive %>% group_by(V2) %>% dplyr::summarise(VerySensitive=n()) %>% merge(res,.,by="V2")
res<-Sensitive %>% group_by(V2) %>% dplyr::summarise(Sensitive=n()) %>% merge(res,.,by="V2")
res<-Specific %>% group_by(V2) %>% dplyr::summarise(Specific=n()) %>% merge(res,.,by="V2")
res<-res %>% mutate(VerySensitive=VerySensitive/Total) %>% mutate(Sensitive=Sensitive/Total) %>% mutate(Specific=Specific/Total)
colnames(res)[1]<-"Cell"
return(res)
}


## Functions for redeem Plus  2024-8-5

add_freq_raw<-function(raw){
GiveName <- c("UMI", "Cell", "Pos", "Variants", "Call", "Ref", "FamSize", "GT_Cts", "CSS", "DB_Cts", "SG_Cts", "Plus", "Minus", "Depth")
colnames(raw)<-GiveName
raw$CellVar<-paste(raw$Cell,raw$Variants,sep="_")
raw.gtsummary<-GTSummary(raw)
raw<-merge(raw,raw.gtsummary[,c("Var1","Freq")],by.x="CellVar",by.y="Var1",all.x = T)
return(raw)
}

#' Function needed to compute the 
make_position_df_3.4<-function(in_df){
    first <- str_split_fixed(in_df[,"UMI"], "_", 3)[, c(2)] %>% 
        as.numeric()
    last <- str_split_fixed(in_df[,"UMI"], "_", 3)[, c(3)] %>% 
        as.numeric()
    start <- pmin(first, last)
    end <- pmax(first, last)
    df <- data.frame(UMI=in_df[,"UMI"],start = start, end = end, pos = in_df$Pos, 
        variant = in_df$Variants) %>% mutate(length = end - start) %>% 
        mutate(rel_position = (pos - start)/length,Freq=in_df[,"Freq"])
    df$edge_dist<-pmin(abs(df$pos-df$start),abs(df$pos-df$end))
    return(df)
}


#' Produce a raw fragment table with frequency (how many cells) and the reletive distance
#' from redeemR object, 
#' 
#' @param redeemR  a redeemR object 
#' @export
add_raw_fragment <- function(redeemR,raw="RawGenotypes.Sensitive.StrandBalance"){
    redeemR.raw <- read.table(paste(redeemR@attr$path,raw,sep="/"))
    filtered.variants <- unique(redeemR@GTsummary.filtered$Variants)
    redeemR.raw.passfilter <- subset(redeemR.raw, V4 %in% filtered.variants)
    raw.pos<- redeemR.raw.passfilter %>% add_freq_raw() %>% make_position_df_3.4()
    return(add_raw_fragment)
}

#' Convinient function that takes raw fragment, and out put fragment with frequency
add_freq_raw<-function(raw){
GiveName <- c("UMI", "Cell", "Pos", "Variants", "Call", "Ref", "FamSize", "GT_Cts", "CSS", "DB_Cts", "SG_Cts", "Plus", "Minus", "Depth")
colnames(raw)<-GiveName
raw$CellVar<-paste(raw$Cell,raw$Variants,sep="_")
raw.gtsummary<-GTSummary(raw)
raw<-merge(raw,raw.gtsummary[,c("Var1","Freq")],by.x="CellVar",by.y="Var1",all.x = T)
return(raw)
}