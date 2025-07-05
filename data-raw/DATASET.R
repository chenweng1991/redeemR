## code to prepare `DATASET` dataset goes here
# Simple reverse complement function
library(data.table)
library(stringi)
library(dplyr)
library(gsubfn)
##Prepare Mito sequence context
mitoref<-read.table("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/Run210706L2/COMBINE_CALLV/CW_mgatk_test/final/chrM_refAllele.txt")
mitobasebias<-table(toupper(mitoref$V2)) %>% as.data.frame()
# Important! Make Context dictionary
# library("insect",lib="/home/cweng/R/x86_64-pc-linux-gnu-library/4.1-focal")


rc<-function(input){
o<-stri_reverse(gsubfn(".", list("A" = "T", "T" = "A","C"="G", "G"="C"), input))
return(o)
}

Contexts<-c()
for(i in 2:(nrow(mitoref)-1)){
    Context<-paste(mitoref$V2[i-1],mitoref$V2[i],mitoref$V2[i+1],sep="")
    Contexts<-c(Contexts,toupper(Context))
}
ContextsDic<-c("GGA",Contexts,"TGG")
names(ContextsDic)<-as.character(1:nrow(mitoref))
head(ContextsDic)
idx<-which(substr(ContextsDic,2,2) %in% c("G","A"))
ContextsDic[idx]<-stri_reverse(rc(ContextsDic[idx]))

# This is actually just complement. It's ok.
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}

ref_all <- fread("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/Run210706L2/COMBINE_CALLV/CW_mgatk_test/final/chrM_refAllele.txt")
colnames(ref_all) <- c("pos", "ref")
ref_all$ref <- toupper(ref_all$ref)
l <- as.character(ref_all$ref)

# Gs happen to be at the first and last position
ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

# Remove Ns
ref_all <- ref_all[!grepl("N", ref_all$three),]

# Make every possible mutation
ref_all_long <- rbind(ref_all,ref_all, ref_all,ref_all)
ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]

# add some meta data
ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

# A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
table(ref_all$ref) # so the reference strand is light (more C/T)
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")

# Change to C/T as ref allele
ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)


usethis::use_data(ref_all_long,ContextsDic, overwrite = TRUE)


# CellPCT
# CellPCT$muRate<-CellPCT[,c("PCT.D4BM","PCT.D4HPC","PCT.D4HSC","PCT.D9BM","PCT.D9HPC","PCT.D9HSC","PCT.D1BM","PCT.D1HPC")] %>% apply(.,1,median)
CellPCT<-readRDS("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9/RDS/DN1_DN4_DN9.CellPCT")
usethis::use_data(CellPCT,overwrite = TRUE)

library(EZsinglecell)
data(msig.db)
data(all.genes.refer)
usethis::use_data(msig.db,overwrite = TRUE)
usethis::use_data(all.genes.refer,overwrite = TRUE)

## 2024-8-11 add more infor into the CellPCT
Old1_BMMC_mitoTracing.sensitive<-readRDS("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9_full/RDS/Old1_BMMC_mitoTracing.sensitive")
Old1_HSPC_mitoTracing.sensitive<-readRDS("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9_full/RDS/Old1_HSPC_mitoTracing.sensitive")
Old2_BMMC_mitoTracing.sensitive<-readRDS("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9_full/RDS/Old2_BMMC_mitoTracing.sensitive")
Old2_HSPC_mitoTracing.sensitive<-readRDS("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9_full/RDS/Old2_HSPC_mitoTracing.sensitive")

Old1_BMMC_redeemR.sensitive<-convert_mitotracing_redeemR(Old1_BMMC_mitoTracing.sensitive,addAssignedVariant = F)
Old1_HSPC_redeemR.sensitive<-convert_mitotracing_redeemR(Old1_HSPC_mitoTracing.sensitive,addAssignedVariant = F)
Old2_BMMC_redeemR.sensitive<-convert_mitotracing_redeemR(Old2_BMMC_mitoTracing.sensitive,addAssignedVariant = F)
Old2_HSPC_redeemR.sensitive<-convert_mitotracing_redeemR(Old2_HSPC_mitoTracing.sensitive,addAssignedVariant = F)

Freq.old1_BM<-Old1_BMMC_redeemR.sensitive@V.fitered[,c("Variants","CellNPCT")] %>% rename(Old1_BM=CellNPCT)
Freq.old1_HSPC<-Old1_HSPC_redeemR.sensitive@V.fitered[,c("Variants","CellNPCT")] %>% rename(Old1_HSPC=CellNPCT)
Freq.old2_BM<-Old2_BMMC_redeemR.sensitive@V.fitered[,c("Variants","CellNPCT")] %>% rename(Old2_BM=CellNPCT)
Freq.old2_HSPC<-Old2_HSPC_redeemR.sensitive@V.fitered[,c("Variants","CellNPCT")] %>% rename(Old2_HSPC=CellNPCT)

CellPCT.update<-merge(CellPCT,Freq.old1_BM, all=T) %>% merge(.,Freq.old1_HSPC, all=T) %>% 
merge(.,Freq.old2_BM, all=T) %>% merge(.,Freq.old2_HSPC, all=T)
CellPCT.update[is.na(CellPCT.update)]<-0

usethis::use_data(CellPCT.update,overwrite = TRUE)


## 2025-07-02 variant annotation
# Mitomap frequency data
library(openxlsx)
mitomap_freq <- read.xlsx("marker_finder.xlsx",
                startRow = 7,
                colNames = TRUE)
mitomap_freq$Overall.Variant.Frequency.in.Sequence.Set   <- as.numeric(sub("%", "", mitomap_freq$Overall.Variant.Frequency.in.Sequence.Set))
mitomap_freq$Frequency.in.Lineage.L      <- as.numeric(sub("%", "", mitomap_freq$Frequency.in.Lineage.L))
mitomap_freq$Frequency.in.Lineage.M      <- as.numeric(sub("%", "", mitomap_freq$Frequency.in.Lineage.M))
mitomap_freq$Frequency.in.Lineage.N      <- as.numeric(sub("%", "", mitomap_freq$Frequency.in.Lineage.N))
mitomap_freq<- mitomap_freq %>% mutate(ID = paste(Position, rCRS.nt, Variant.nt,sep = "_"))
names(mitomap_freq)[12] <- "RSRS50"

# Mitomap haplogroup markers
haplogroup_markers <- read.table("freq_variants.txt", header = T)
haplogroup_markers <- haplogroup_markers %>% mutate(ID = paste(tpos, tnt, qnt, sep = "_"))

# mitochondrial diseases
mito_diseases <- read.csv("MutationsCodingControl_MITOMAP.csv", header = T)
mito_diseases$ID <- paste(mito_diseases$Position,sub("-","_", mito_diseases$NucleotideChange), sep = "_")

# homopolymer bed
mito_homopolymer <- read.table("/lab/solexa_weissman/cweng/Projects/Collaborator/Sankaran_lineage_tracing/ReDeeM2.0/01_mitoblacklist_script/mt_homopolymers_4bp.bed", header = F)

usethis::use_data(mitomap_freq,overwrite = TRUE)
usethis::use_data(haplogroup_markers,overwrite = TRUE)
usethis::use_data(mito_diseases,overwrite = TRUE)
usethis::use_data(mito_homopolymer,overwrite = TRUE)

# 2025-07-04 
file.copy(
  from      = "/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Revision20230904/Reviewer5_data/GRCH38_MT_bioMRT_output_refcds.rda",
  to        = "../data/GRCH38_MT_bioMRT_output_refcds.rda",
  overwrite = TRUE    # overwrite if it already exists
)
