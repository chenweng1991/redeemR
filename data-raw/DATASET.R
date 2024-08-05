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