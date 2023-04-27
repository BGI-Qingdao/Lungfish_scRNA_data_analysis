library(Seurat)
library(dplyr)
library(cowplot)
library(reshape2)
library(DoubletFinder)
library(data.table)

args<-commandArgs(T)
if(length(args)!=5){
print ("Usage:")
print ("cd work_dir ;Rcript anno.R count_matrix sample max_nFeature_RNA doublet_rate")
print ("Note:")
print ("Aim to run single-sample umap")
q(save = "no", status = 0, runLast = TRUE)
}
setwd(args[1])

### 1.Create object
mtx <- Read10X(args[2], gene.column = 1)
object.data = mtx
sample<-args[3]
colnames(object.data)<-paste(colnames(object.data),sample,sep="_")
object_name <- CreateSeuratObject(counts = object.data, project =sample, min.cells = 3, min.features = 200)

if (sum(object_name@meta.data$percent.mt)== 0){
            object_name@meta.data$percent.mt[1]<-0.000001
}
pdf(paste0(sample,"raw.vln.pdf"),width=15,height=10)
VlnPlot(object_name,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()
object_name <- subset(x = object_name, subset = nFeature_RNA > 200 & nFeature_RNA < as.numeric(args[4]))
a=length(colnames(object_name))
print(a)
pdf(paste0(sample,"filterGene.vln.pdf"),width=15,height=10)
VlnPlot(object_name,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()

### 2.normalize
object_name<- SCTransform(object_name, verbose = FALSE,vars.to.regress="percent.mt")
object_name <- RunPCA(object_name, features = VariableFeatures(object = object_name))
object_name <- FindNeighbors(object_name, dims = 1:30)
object_name <- RunUMAP(object_name, dims = 1:30)
object_name <- FindClusters(object_name, resolution = 0.5)

### 3.doublet
sweep.res.list_SM <- paramSweep_v3(object_name, PCs = 1:30,sct=TRUE)
sweep.stats_SM <- summarizeSweep(sweep.res.list_SM)
bcmvn_SM <- find.pK(sweep.stats_SM)
pK_value <- as.numeric(as.character(bcmvn_SM$pK[bcmvn_SM$BCmetric == max(bcmvn_SM$BCmetric)]))
annotations <- object_name@meta.data$RNA_snn_res.1
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(as.numeric(args[5])*length(object_name@meta.data$orig.ident))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
pN_value <- 0.25
pANN_value <- paste0("pANN_",pN_value,"_",pK_value,'_',nExp_poi)
object_name <- doubletFinder_v3(object_name, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi, reuse.pANN = FALSE,sct=TRUE)
object_name <- doubletFinder_v3(object_name, PCs = 1:30, pN = pN_value, pK = pK_value, nExp = nExp_poi.adj, reuse.pANN = pANN_value,sct=TRUE)
object_name@meta.data[,"DF_hi.lo"] <- object_name@meta.data[,10]
object_name@meta.data$DF_hi.lo[which(object_name@meta.data$DF_hi.lo == "Doublet" & object_name@meta.data[,11] == "Singlet")] <- "Doublet_lo"
object_name@meta.data$DF_hi.lo[which(object_name@meta.data$DF_hi.lo == "Doublet")] <- "Doublet_hi"
object_name@meta.data$Doublet <-eval(parse(text = paste0("object_name@meta.data$DF.classifications_",pN_value,"_",pK_value,'_',nExp_poi))) 

### 4.subset&plot
object_Singlet <- SubsetData(object_name,subset.name='Doublet',accept.value='Singlet')
object_Singlet<-SCTransform(object_Singlet, verbose = FALSE,vars.to.regress="percent.mt")
pdf("Vlnplot.pdf",width=15,height=10)
VlnPlot(object_Singlet,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()
object_Singlet <- RunPCA(object = object_Singlet, npcs = , verbose = FALSE)
pdf("DimHeatmap_30pc.pdf",width=15,height=10)
DimHeatmap(object = object_Singlet, dims = 1:30, cells= 200, balanced = TRUE)
dev.off()
object_Singlet <- RunUMAP(object = object_Singlet, reduction = "pca", dims = 1:30)
object_Singlet <- FindNeighbors(object = object_Singlet, reduction = "pca", dims = 1:30)
object_Singlet <- FindClusters(object_Singlet, resolution = 0.5)
pdf("sample-umap.pdf",width=15,height=10)
DimPlot(object =object_Singlet, reduction = "umap")
dev.off()
object_Singlet=NormalizeData(object = object_Singlet, assay ="RNA",normalization.method = "LogNormalize",scale.factor = 10000)
object_Singlet=ScaleData(object_Singlet, verbose = FALSE,assay="RNA")
saveRDS(object_Singlet, file = "combined_analysis.rds")
