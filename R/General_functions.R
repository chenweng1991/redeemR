#' Tomerge_v2
#'
#' This function is to quickly merge two dataframe by rownames, but can choose to leave A or B all information
#' @param A  dataframe A
#' @param B  dataframe B
#' @return  return a data frame with merged information
#' @export
#' @examples
#' Tomerge_v2(A,B)

Tomerge_v2<-function(A,B,leavex=T,leavey=F){
	mergeAB<-merge(A,B,by="row.names",all.x=leavex,all.y=leavey)
	row.names(mergeAB)<-mergeAB[,1]
	mergeAB<-mergeAB[,-1]
	return(mergeAB)
}



#' Function to translate the RNA barcode into ATAC barcode and add a column
#'
#' This function allows you to input the metadata with row name as cell barcode
#' @param meta  a dataframe with the row names as the RNA cell barcode usually with the post -1
#' @param bclength The cell barcode length, default is 16
#' @param from A vector of the postfix,  usually is c(1,2,3,...), it depends on how many samples are aggregated in Cellranger RNA part
#' @param to A vector of the postfix, those cooresponds to the postfix added in redeemR, in general, if it matches, then simply c(1,2,3,...), 
#' but in case not match, here provides a way to transform into scredeemR order
#' @return meta a dataframe
#' @examples
#' Translate_RNA2ATAC(meta)
#' @export
Translate_RNA2ATAC<-function(meta=bmmc.filtered@meta.data,PostFix=T,bclength=16,from=c(1,2,3),to=c(1,2,3)){
data(ATACWhite)
data(RNAWhite)
L<-nchar(row.names(meta))
post<-substr(row.names(meta),bclength+2,L)
post<-post %>% plyr::mapvalues(.,from=from,to=to)
NakeName<-substr(row.names(meta),1,bclength)
Dic2<-ATACWhite$V1
names(Dic2)<-as.character(RNAWhite$V1)
if(PostFix){
	meta$ATACName<-paste(Dic2[NakeName],post,sep="_")
}else{
	meta$ATACName<-Dic2[NakeName]
}
return(meta)
}

#' Translate_simple_RNA2ATAC 
#'
#' This function allows you to input the RNA name to translate to ATAC name 
#' @param name  RNA name, as the RNA cell barcode usually with the post -1
#' @param bclength The cell barcode length, default is 16
#' @param from A vector of the postfix,  usually is c(1,2,3,...), it depends on how many samples are aggregated in Cellranger RNA part
#' @param to A vector of the postfix, those cooresponds to the postfix added in redeemR, in general, if it matches, then simply c(1,2,3,...), 
#' but in case not match, here provides a way to transform into redeemR order
#' @return ATAC name
#' Translate_RNA2ATAC(`a vector of RNA names`)
#' @export
Translate_simple_RNA2ATAC<-function(name,PostFix=T,bclength=16,from=c(1,2,3),to=c(1,2,3)){
data(ATACWhite)
data(RNAWhite)
L<-nchar(name)
post<-substr(name,bclength+2,L)
post<-post %>% plyr::mapvalues(.,from=from,to=to)
NakeName<-substr(name,1,bclength)
Dic2<-ATACWhite$V1
names(Dic2)<-as.character(RNAWhite$V1)
if(PostFix){
	name<-paste(Dic2[NakeName],post,sep="_")
}else{
	name<-Dic2[NakeName]
}
return(name)
}

#' Translate_simple_ATAC2RNA 
#'
#' This function allows you to input the ATAC name to translate to RNA name 
#' @param name  RNA name, as the RNA cell barcode usually with the post -1
#' @param bclength The cell barcode length, default is 16
#' @param from A vector of the postfix,  usually is c(1,2,3,...), it depends on how many samples are aggregated in Cellranger RNA part
#' @param to A vector of the postfix, those cooresponds to the postfix added in redeemR, in general, if it matches, then simply c(1,2,3,...), 
#' but in case not match, here provides a way to transform into redeemR order
#' @return RNA name
#' Translate_RNA2ATAC(`a vector of RNA names`)
#' @export
Translate_simple_ATAC2RNA<-function(name,PostFix=T,bclength=16,from=c(1,2,3),to=c(1,2,3)){
data(ATACWhite)
data(RNAWhite)
L<-nchar(name)
post<-substr(name,bclength+2,L)
post<-post %>% plyr::mapvalues(.,from=from,to=to)
NakeName<-substr(name,1,bclength)
Dic2<-RNAWhite$V1
names(Dic2)<-as.character(ATACWhite$V1)
if(PostFix){
	name<-paste(Dic2[NakeName],post,sep="-")
}else{
	name<-Dic2[NakeName]
}
return(name)
}





#' Function to Merge sparse Matrix
#'
#' This function allows you to input a list of sparse matrix and merge by rownames, return a new sparse matrix
#' @param mtx.list  A list of sparse matrix to be merged
#' @param postfix  a vector of postfix (Usually are numbers that added at the end of cell names). Better be consistent with a merged redeemR object orders
#' @return new sparse matrix
#' @examples
#' Donor4_HSC_HPC_BMMC.Mtx<-MergeMtx(list(Donor04_BMMC_Multiome_wrapper$seurat@assays$RNA@counts,Donor04_HPC_Multiome_wrapper$seurat@assays$RNA@counts,Donor04_HSC_Multiome_wrapper$seurat@assays$RNA@counts),c(3,2,1))
#' Donor4_HSC_HPC_BMMC.RNA.seurat<-GEM_Wrapper(Donor4_HSC_HPC_BMMC.Mtx)
#' @export
MergeMtx<-function(mtx.list,postfix){
colnames(mtx.list[[1]])<-strsplit(colnames(mtx.list[[1]]),"-") %>% lapply(.,function(x){x[1]}) %>% unlist %>% paste(.,postfix[1],sep="-")
Merged.mtx<-as.matrix(mtx.list[[1]])
for(i in 2:length(mtx.list)){
    colnames(mtx.list[[i]])<-strsplit(colnames(mtx.list[[i]]),"-") %>% lapply(.,function(x){x[1]}) %>% unlist %>% paste(.,postfix[i],sep="-")
    Merged.mtx<-Tomerge_v2(Merged.mtx,as.matrix(mtx.list[[i]]),leavex = T, leavey = T)
}
Merged.mtx[is.na(Merged.mtx)]<-0
Merged.mtx<-Matrix(as.matrix(Merged.mtx))
return(Merged.mtx)
}

#' Function to add hematopoietic signatures from Griffin_Signatures
#'
#' This function allows you to input a seurat object, add the signatures and return an seurat object
#' @param object a seurat object
#' @return a seurat object
#' @export
AddHemSignature<-function(object=Donor01_BMMC_Multiome_wrapper.filtered){
require(Seurat)
data(Griffin_Signatures)
Sig.HSC<-list(Griffin_Signatures$Griffin_BPDCN_HSC %>% .[!is.na(.)])
Sig.Prog<-list(Griffin_Signatures$Griffin_BPDCN_Prog %>% .[!is.na(.)])
Sig.EarlyE<-list(Griffin_Signatures$Griffin_BPDCN_EarlyE %>% .[!is.na(.)])
Sig.LateE<-list(Griffin_Signatures$Griffin_BPDCN_LateE	 %>% .[!is.na(.)])
Sig.ProMono<-list(Griffin_Signatures$Griffin_BPDCN_ProMono %>% .[!is.na(.)])
Sig.Mono<-list(Griffin_Signatures$Griffin_BPDCN_Mono %>% .[!is.na(.)])
Sig.ncMono<-list(Griffin_Signatures$Griffin_BPDCN_ncMono %>% .[!is.na(.)])
Sig.cDC<-list(Griffin_Signatures$Griffin_BPDCN_cDC %>% .[!is.na(.)])
Sig.pDC<-list(Griffin_Signatures$Griffin_BPDCN_pDC %>% .[!is.na(.)])
Sig.ProB<-list(Griffin_Signatures$Griffin_BPDCN_ProB %>% .[!is.na(.)])
Sig.PreB<-list(Griffin_Signatures$Griffin_BPDCN_PreB %>% .[!is.na(.)])
Sig.B<-list(Griffin_Signatures$Griffin_BPDCN_B %>% .[!is.na(.)])
Sig.Plasma<-list(Griffin_Signatures$Griffin_BPDCN_Plasma %>% .[!is.na(.)])
Sig.T<-list(Griffin_Signatures$Griffin_BPDCN_T %>% .[!is.na(.)])
Sig.CTL<-list(Griffin_Signatures$Griffin_BPDCN_CTL %>% .[!is.na(.)])
Sig.NK<-list(Griffin_Signatures$Griffin_BPDCN_NK %>% .[!is.na(.)])
DefaultAssay(object)<-"SCT"
object<-AddModuleScore(object = object, features = Sig.HSC, name = "Sig.HSC")
object<-AddModuleScore(object = object, features = Sig.Prog, name = "Sig.Prog")
object<-AddModuleScore(object = object, features = Sig.EarlyE, name = "Sig.EarlyE")
object<-AddModuleScore(object = object, features = Sig.LateE, name = "Sig.LateE")
object<-AddModuleScore(object = object, features = Sig.ProMono, name = "Sig.ProMono")
object<-AddModuleScore(object = object, features = Sig.Mono, name = "Sig.Mono")
object<-AddModuleScore(object = object, features = Sig.ncMono, name = "Sig.ncMono")
object<-AddModuleScore(object = object, features = Sig.cDC, name = "Sig.cDC")
object<-AddModuleScore(object = object, features = Sig.pDC, name = "Sig.pDC")
object<-AddModuleScore(object = object, features = Sig.ProB, name = "Sig.ProB")
object<-AddModuleScore(object = object, features = Sig.PreB, name = "Sig.PreB")
object<-AddModuleScore(object = object, features = Sig.B, name = "Sig.B")
object<-AddModuleScore(object = object, features = Sig.Plasma, name = "Sig.Plasma")
object<-AddModuleScore(object = object, features = Sig.T, name = "Sig.T")
object<-AddModuleScore(object = object, features = Sig.CTL, name = "Sig.CTL")
object<-AddModuleScore(object = object, features = Sig.NK, name = "Sig.NK")
return(object)    
}


#' Function to reclustering a seurat object
#'
#' This function allows you to input a seurat object(multiome), redo clustering. Usually this is after subset
#' @param ob a seurat object
#' @return a seurat object
#' @export
Reclustering<-function(ob){
DefaultAssay(ob) <- "RNA"
ob <- SCTransform(ob, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
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

#' Function to reclustering_hm a seurat object with Harmony
#'
#' This function allows you to input a seurat object(multiome), redo clustering harmony by a certain column in meta data. Usually this is after subset
#' @param ob a seurat object
#' @param HarmonyBy The columne name in meta that will be used for Harmony
#' @return a seurat object
#' @export
Reclustering_hm<-function(ob=DN4_RigHSC_T1T2_Multiome_wrapper_filtered.anno,HarmonyBy="TimePoint"){
DefaultAssay(ob) <- "RNA"
ob <- suppressWarnings(SCTransform(ob, verbose = FALSE) %>% RunPCA() )
ob<- RunHarmony(ob, group.by.vars = HarmonyBy,assay.use = "SCT",reduction.save = "harmony.rna",project.dim = F)
ob<-RunUMAP(ob,dims = 1:50, reduction.name = 'umap.H.rna',reduction = "harmony.rna", reduction.key = 'HrnaUMAP_')
DefaultAssay(ob) <- "ATAC"
ob <- RunTFIDF(ob)
ob <- FindTopFeatures(ob, min.cutoff = 'q0')
ob <- RunSVD(ob)
ob<- RunHarmony(ob, group.by.vars = "TimePoint",assay.use = "ATAC",reduction.save = "harmony.atac",project.dim = F)
ob <- RunUMAP(ob, reduction = 'harmony.atac', dims = 2:50, reduction.name = "umap.H.atac", reduction.key = "HatacUMAP_")
ob <- FindMultiModalNeighbors(ob, reduction.list = list("harmony.rna", "harmony.atac"), dims.list = list(1:50, 2:50))
ob <- RunUMAP(ob, nn.name = "weighted.nn", reduction.name = "umap.H.wnn", reduction.key = "HwnnUMAP_")
ob <- FindClusters(ob, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
return(ob)
}


##----------------------------------------------Below are two DE functions from EZsinglecell-------------------------------------------------------##
#' DE.gettripple
#'
#' This function is to prepare the data format that is used to differentially expression calling. It include the raw matrix; data.info and size effect
#' @param datapair    tyhe datapair generated from datapair.mk
#' @param cpcol The column name for comparison.
#' @param withscran  if true, use deconvolution to calculate size effect.
#' @return  This will return .tri.dummy file that is the input for DE analysis
#' @export
#' @examples
#' ROCKvsnorock.endo.tri.dummy<-DE.gettripple(ROCKvsnorock.endo.paired,cpcol="name")
DE.gettripple<-function(datapair,cpcol,withscran=F)
{
#	library(scran)
	x<-datapair$data
	x[is.na(x)]<-0
	x.info<-datapair$info
	if(!all(colnames(x) ==row.names(x.info)))
	{
		stop("The datapair is not match well")
	}else
	{
		print("datapair match well......continue")
	}
#
	wherezero<-function(vector)
	{

		if (all(vector==0))
		{
			return(FALSE)
		}
		else
		{
			return(TRUE)
		}
	}
#
	makedummycell<-function(vector)
	{
		dummycell<-0
		if(all(vector==0))
		{
			dummycell<-1
		}
		return(dummycell)
	}
#
	x<-x[,apply(x,2,wherezero)]  # get rid of cells with no expression
	if(withscran)
	{
		sf<-computeSumFactors(as.matrix(x), positive=F)  # calculate  size factor
		x<-x[,which(sf>0)]
		x.info<-x.info[which(sf>0),]
		sf<-sf[which(sf>0)]
	}else
	{
		sf<-colSums(x)
	}
	for (celltype in unique(x.info[,cpcol]) )
	{
		subdata<-x[,row.names(x.info)[which(x.info[,cpcol]==celltype)]]
		subdata[is.na(subdata)]<-0
		x<-cbind(x,apply(subdata,1,makedummycell))
		colnames(x)[ncol(x)]<-celltype
		newinfoline<-x.info[1,]
		newinfoline[,cpcol]<-celltype
		x.info<-rbind(x.info,newinfoline)
		row.names(x.info)[nrow(x.info)]<-celltype
		sf<-c(sf,median(sf))
		names(sf)[length(sf)]<-celltype
	}
	return(list(data=x,info=x.info,sf=sf))
}


#' DoDE
#'
#' This is the main function for calculating differentially expressed genes
#' @param tri.dummy this is generated from  DE.gettripple
#' @param cpcol  the column in tri.dummy$info, the contents of which are used for iteratively compare with one another
#' @param onlyoneSample If true, regress out batch effect. Notice, there should be a "Sample" column in in tri.dummy$info that indicate sample or donor or batch
#' @param cpus a number of cpus being used for calculation, default is 16
#' @return  return a list that includes all DE result iteratively
#' @export
#' @examples
#' ROCKvsnorock.endo.de<-DoDE(ROCKvsnorock.endo.tri.dummy,"name",onlyoneSample=T,cpus=16)

DoDE<-function(tri.dummy,cpcol,onlyoneSample=F,cpus=16)
{

Donbregression<-function(data,info,cpcol,reftype,sf,gene,onlyoneSample){
		require(MASS)
		info[,cpcol]<-factor(info[,cpcol],levels=c(reftype,setdiff(levels(info[,cpcol]),reftype)))
		if(onlyoneSample)
		{
			md<-try(glm.nb(data[gene,] ~ info[,cpcol]+log(sf)),silent=T)
		}else
		{
			md<-try(glm.nb(data[gene,] ~ info[,cpcol]+info[,"Sample"]+log(sf)),silent=T)
		}
					coeff.table<-try(summary(md)$coefficients,silent=T)
					N<-2
					output<-c()
					for(i in setdiff(levels(info[,cpcol]),reftype))
					{
							p.value<-try(coeff.table[N,4],silent=T)
							if(grepl("Error",p.value))
							{
									p.value<-NA
							}
				ref.mean<-mean(data[,row.names(info)[which(info[,cpcol]==reftype)]])
				i.mean<-mean(data[,row.names(info)[which(info[,cpcol]==i)]])
							assign(gene,data.frame(
							ref.sum=as.integer(sum(data[gene,row.names(info)[which(info[,cpcol]==reftype)]])),
							alt.sum=as.integer(sum(data[gene,row.names(info)[which(info[,cpcol]==i)]])),
							ref.nonzero=as.integer(length(which(data[gene,row.names(info)[which(info[,cpcol]==reftype)]]>0))),
							alt.nonzero=as.integer(length(which(data[gene,row.names(info)[which(info[,cpcol]==i)]]>0))),
							nonzero.ratio=(length(which(data[gene,row.names(info)[which(info[,cpcol]==reftype)]]>0))/as.integer(as.matrix(table(info[,cpcol]))[reftype,]))/(length(which(data[gene,row.names(info)[which(info   [,cpcol]==i)]]>0))/as.integer(as.matrix(table(info[,cpcol]))[i,])),
							mean.ratio=(mean(data[gene,row.names(info)[which(info[,cpcol]==reftype)]])/ref.mean)/(mean(data[gene,row.names(info)[which(info[,cpcol]==i)]])/i.mean),
							p.value=p.value
							))
							cline<-get(gene)
							row.names(cline)<-gene
							clist<-list(cline)
							names(clist)<-i
							output<-c(output,clist)
							N<-N+1
					}
					return(output)
	}

calculateDE.sg<-function(data,info,cpcol,sf,gene,onlyoneSample){
	singlegene.result<-list()
	length(singlegene.result)<-length(levels(info[,cpcol]))
	names(singlegene.result)<-levels(info[,cpcol])
	for (reftype in levels(info[,cpcol])){
		singlegene.result[[reftype]]<-Donbregression(data,info,cpcol,reftype,sf,gene,onlyoneSample)
		}
		return(singlegene.result)
}

rbindlist<-function(pardoresult,info,cpcol){
		print("Start binding......")
		DEresult<-list()
		length(DEresult)<-length(levels(info[,cpcol]))
		names(DEresult)<-levels(info[,cpcol])
		for (reftype in levels(info[,cpcol]))
		{
			DEresult[[reftype]]<-list(DEresult[[reftype]])
			length(DEresult[[reftype]])<-length(levels(info[,cpcol]))-1
			names(DEresult[[reftype]])<-setdiff(levels(info[,cpcol]),reftype)
		}
		n<-0
		for(reftype in levels(info[,cpcol]))
		{
			print(n)
			n<-n+1
	    		for (alttype in setdiff(levels(info[,cpcol]),reftype))
	    		{
	        		for (i in 1:length(pardoresult))
	        		{
	            			DEresult[[reftype]][[alttype]]<-rbind(DEresult[[reftype]][[alttype]],pardoresult[[i]][[reftype]][[alttype]])
	        		}
	    		}
		cellnumber<-list(as.matrix(table(info[,cpcol]))[reftype,])
		names(cellnumber)="cellnumber"
		DEresult[[reftype]]<-c(DEresult[[reftype]],cellnumber)
		}
		return(DEresult)
	}

# main function starting from here
data<-tri.dummy[[1]]
info<-tri.dummy[[2]]
sf<-tri.dummy[[3]]
data<-as.matrix(data)
if(!all(colnames(data) ==row.names(info)))
{
	stop("datapair doesnt match well....data:info doesn't match")
}else
{
	print("data::info match well......continue")
}
if(!all(colnames(data) ==row.names(sf)))
{
	print("datapair doesnt match well....data:sf doesn't match")
}else
{
	print("data::sf match well......continue")
}
print(paste("The compared column is",cpcol))
info[,cpcol]<-factor(info[,cpcol],levels=as.character(unique(info[,cpcol])))
print(table(info[,cpcol]))
timestart<-Sys.time()
require(doParallel)
cl <- makeCluster(cpus)
doParallel::registerDoParallel(cl)
time.start<-Sys.time()
for (reftype in levels(info[,cpcol]))
{
	assign(reftype,list())
}
pardoresult<-list()
pardoresult<-foreach (gene=row.names(data)) %dopar%
{
	print(which(row.names(data)==gene))
	print(calculateDE.sg(data,info,cpcol,sf,gene,onlyoneSample))
}
print("pardoresult Get!!")
time.ends<-Sys.time()
DEresult<-rbindlist(pardoresult,info,cpcol)
timeend<-Sys.time()
print(timeend-timestart)
return(DEresult)
}

#' Motifenrich.binom 
#' In house function to compute enrichment from Fimo
#' This function was developed based on HSC_multiome_Het.ipynb and HSC_multiome_Het_2.ipynb
#' @param queryP.motif        can be a subset of all.motif.sig
#' @param controlP.motif      can be all.motif.sig
#' @param alt default is greater
#' @export
#' @import qvalue
Motifenrich.binom<-function(queryP.motif,controlP.motif,alt="greater"){
require(qvalue)
query.peaks.motif.summary<-table(queryP.motif$motif_id) %>% as.data.frame %>% .[order(.$Freq,decreasing=T),]
all.motif.sig.summary<-table(controlP.motif$motif_id) %>% as.data.frame %>% .[order(.$Freq,decreasing=T),]
query.peaks.motif.summary<-merge(all.motif.sig.summary,query.peaks.motif.summary,by="Var1")
query.peaks.motif.summary<-subset(query.peaks.motif.summary,Freq.x>0)
query.peaks.motif.summary$Freq.x<-query.peaks.motif.summary$Freq.x+1
query.peaks.motif.summary$Freq.y<-query.peaks.motif.summary$Freq.y+1
query.peaks.motif.summary<- cbind(query.peaks.motif.summary,ratio=query.peaks.motif.summary[,"Freq.y"]/query.peaks.motif.summary[,"Freq.x"])
pvalues<-c()
for(i in 1:nrow(query.peaks.motif.summary)){
	 # print(i)
binom.md<-binom.test(query.peaks.motif.summary[i,3],sum(query.peaks.motif.summary[,3]),query.peaks.motif.summary[i,2]/sum(query.peaks.motif.summary[,2]),alternative=alt)
pvalues<-c(pvalues,binom.md$p.value)
}
fulllist<-cbind(query.peaks.motif.summary,FC=query.peaks.motif.summary$ratio/(sum(query.peaks.motif.summary$Freq.y)/sum(query.peaks.motif.summary$Freq.x)),pvalues=pvalues)
fulllist$qvalue<-qvalue(fulllist$pvalues)$qvalues
return(fulllist)
}

#' Fun.enrich_withFC
#'
#' This function is to do GSEA enrichment analysis
#' @param markergenesList is a list containing several elements each of which is a vector of gene names
#' @param All.genes the background genes  data(all.genes.refer)
#' @param db  =msig.db   the term database to use, default is the msig.db;  data(msig.db)
#' @param qcutoff=0.05  the term to show with q value lower than 0.05
#' @param top=NULL  pass a number to only show top # terms
#' @return this will return a simple object with two elements, one is the table to show the significant terms , the other one shows the gene names involved significant terms
#' @export
#' @examples
#' S7G.GSEA<-Fun.enrich_withFC(markergenesList=list(S7G=topS7G.genes),All.genes=all.genes.refer)
Fun.enrich_withFC<-function(markergenesList,All.genes=all.genes.refer
,db=msig.db,qcutoff=0.05,top=NULL)
{
	require(qvalue)
	binomial.test.resultsList<-list()
	resultname<-c()
	for (markname in names(markergenesList))
	{
		if(length((markergenesList[[markname]]))==0)
		{
			next
		}
		print(markname)
		binomialtest.msig.result<-c()
		for (name in names(db))
		{
			binomialtest.msig.result<-rbind(binomialtest.msig.result,binomialtest.msig.enrch_deplet(as.character(markergenesList[[markname]]),All.genes,name,thedatabase=db))
		}
		fdr.enrich<-qvalue(binomialtest.msig.result[,2])$qvalues
		fdr.deplete<-qvalue(binomialtest.msig.result[,3])$qvalues
		FDR.combined<-c()
		FDR.combined[which(is.na(binomialtest.msig.result$ratio))]<-NA
		FDR.combined[which(binomialtest.msig.result$ratio>1)]<-fdr.enrich[which(binomialtest.msig.result$ratio>1)]
		FDR.combined[which(binomialtest.msig.result$ratio<=1)]<-fdr.deplete[which(binomialtest.msig.result$ratio<=1)]
		if(is.null(FDR.combined))
		{
			next
		}
		binomialtest.msig.result<-cbind(binomialtest.msig.result,FDR=FDR.combined)
		binomialtest.msig.result<-binomialtest.msig.result[order(binomialtest.msig.result$FDR),]
		binomial.test.resultsList<-c(binomial.test.resultsList,list(binomialtest.msig.result))
		resultname<-c(resultname,markname)

	}

	names(binomial.test.resultsList)=resultname
	DRscombo<-data.frame()
	for (markname in names(binomial.test.resultsList))
	{
		tmp<-subset(binomial.test.resultsList[[markname]],FDR<qcutoff)
		row.names(tmp)<-tmp[,"name"]
		colnames(tmp)[c(4,5)]<-paste(markname,colnames(tmp)[c(4,5)],sep="_")
		for (othermarkname in setdiff(names(binomial.test.resultsList),markname))
		{
			tmp2<-binomial.test.resultsList[[othermarkname]]
			row.names(tmp2)<-tmp2[,"name"]
			tmp2<-tmp2[row.names(tmp),c("ratio","FDR")]
			colnames(tmp2)<-paste(othermarkname,colnames(tmp2),sep="_")
			tmp<-cbind(tmp,tmp2)
		}
		tmp<-tmp[,c(-1,-2,-3)]
		if (!is.null(top))
		{
			tmp<-head(tmp,n=top)
		}
		DRscombo<-Tomerge.col(DRscombo,tmp)
	}
	all.list<-c()
	for (term in row.names(DRscombo))
	{
		cur.genelist<-c()
		for (i in 1:length(markergenesList))
		{
			cur.genes<-markergenesList[[i]][markergenesList[[i]] %in% as.character(db[[term]])]
			cur.genelist<-c(cur.genelist,list(cur.genes))
		}
		names(cur.genelist)<-names(markergenesList)
		all.list<-c(all.list,list(cur.genelist))

	}
	names(all.list)<-row.names(DRscombo)
#)
#DRscombo<-DRscombo[tmpcondition,]
	return(list(DRscombo,all.list=all.list))
}


#' Fun.enrich_withFC.pvalue
#'
#' This function is to do GSEA enrichment analysis, to be updated
#' @param markergenesList is a list containing several elements each of which is a vector of gene names
#' @param All.genes the background genes  data(all.genes.refer)
#' @param db  =msig.db   the term database to use, default is the msig.db;  data(msig.db)
#' @param qcutoff=0.05  the term to show with q value lower than 0.05
#' @param top=NULL  pass a number to only show top # terms
#' @return this will return a simple object with two elements, one is the table to show the significant terms , the other one shows the gene names involved significant terms
#' @export
#' @examples
#'
Fun.enrich_withFC.pvalue<-function(markergenesList,All.genes=all.genes.refer,db=msig.db,qcutoff=0.05,top=NULL)
{
	require(qvalue)
	binomial.test.resultsList<-list()
	resultname<-c()
	for (markname in names(markergenesList))
	{
		if(length((markergenesList[[markname]]))==0)
		{
			next
		}
		print(markname)
		binomialtest.msig.result<-c()
		for (name in names(db))
		{
			binomialtest.msig.result<-rbind(binomialtest.msig.result,binomialtest.msig.enrch_deplet(as.character(markergenesList[[markname]]),All.genes,name,thedatabase=db))
		}
		fdr.enrich<-binomialtest.msig.result[,2]
		fdr.deplete<-binomialtest.msig.result[,3]
		FDR.combined<-c()
		FDR.combined[which(is.na(binomialtest.msig.result$ratio))]<-NA
		FDR.combined[which(binomialtest.msig.result$ratio>1)]<-fdr.enrich[which(binomialtest.msig.result$ratio>1)]
		FDR.combined[which(binomialtest.msig.result$ratio<=1)]<-fdr.deplete[which(binomialtest.msig.result$ratio<=1)]
		if(is.null(FDR.combined))
		{
			next
		}
		binomialtest.msig.result<-cbind(binomialtest.msig.result,FDR=FDR.combined)
		binomialtest.msig.result<-binomialtest.msig.result[order(binomialtest.msig.result$FDR),]
		binomial.test.resultsList<-c(binomial.test.resultsList,list(binomialtest.msig.result))
		resultname<-c(resultname,markname)

	}

	names(binomial.test.resultsList)=resultname
	DRscombo<-data.frame()
	for (markname in names(binomial.test.resultsList))
	{
		tmp<-subset(binomial.test.resultsList[[markname]],FDR<qcutoff)
		row.names(tmp)<-tmp[,"name"]
		colnames(tmp)[c(4,5)]<-paste(markname,colnames(tmp)[c(4,5)],sep="_")
		for (othermarkname in setdiff(names(binomial.test.resultsList),markname))
		{
			tmp2<-binomial.test.resultsList[[othermarkname]]
			row.names(tmp2)<-tmp2[,"name"]
			tmp2<-tmp2[row.names(tmp),c("ratio","FDR")]
			colnames(tmp2)<-paste(othermarkname,colnames(tmp2),sep="_")
			tmp<-cbind(tmp,tmp2)
		}
		tmp<-tmp[,c(-1,-2,-3)]
		if (!is.null(top))
		{
			tmp<-head(tmp,n=top)
		}
		DRscombo<-Tomerge.col(DRscombo,tmp)
	}
	all.list<-c()
	for (term in row.names(DRscombo))
	{
		cur.genelist<-c()
		for (i in 1:length(markergenesList))
		{
			cur.genes<-markergenesList[[i]][markergenesList[[i]] %in% as.character(db[[term]])]
			cur.genelist<-c(cur.genelist,list(cur.genes))
		}
		names(cur.genelist)<-names(markergenesList)
		all.list<-c(all.list,list(cur.genelist))

	}
	names(all.list)<-row.names(DRscombo)
#)
#DRscombo<-DRscombo[tmpcondition,]
	return(list(DRscombo,all.list=all.list))
}



#' binomialtest.msig.enrch_deplet
#'
#' This function is an internal function calculating the significance
#' @param mylist
#' @param All=All.genes
#' @param name
#' @param thedatabase=db
#' @return
#' @export
#' @examples
binomialtest.msig.enrch_deplet<-function(mylist,All=All.genes,name,thedatabase=db)
{
n<-length(mylist)
x<-length(which(mylist %in% thedatabase[[name]]))
p<-length(setdiff(thedatabase[[name]],setdiff(thedatabase[[name]],All)))/length(All)

	binomitest.enrich<-binom.test(x,n,p,alternative="greater")
	binomitest.deplete<-binom.test(x,n,p,alternative="less")
statistics<-data.frame(name,enrichsc=binomitest.enrich$p.value,depletesc=binomitest.deplete$p.value,ratio=(x/n)/p)
return(statistics)
}



#' Tomerge.col
#'
#' This function is an internal function doing ther merge
#' @param df1,
#' @param df2
#' @return ""
#' @export
#' @examples
Tomerge.col<-function(df1,df2)
{
	if(length(df1)==0)
	{
		newdf<-df2
		return(newdf)
	}else
	{
		df2<-df2[!row.names(df2) %in%  row.names(df1),]
		newdf<-data.frame(row.names=c(row.names(df1),row.names(df2)))
		for(name in names(df1))
		{
			if(is.factor(df1[,name]))
			{
				assign(name,c(as.character(df1[,name]),as.character(df2[,name])))
			}else
			{
				assign(name,c(df1[,name],df2[,name]))
			}
			newdf<-data.frame(newdf,get(name))
			names(newdf)[ncol(newdf)]<-name
		}
		return(newdf)
	}
}