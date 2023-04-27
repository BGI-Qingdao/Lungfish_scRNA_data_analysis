library(Seurat)
library(dplyr)
library(cowplot)
library(reshape2)
library(future)
library(harmony)
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 900000 * 1024^2)
args = commandArgs(T)
dir.create(args[1])
setwd(args[1])
method = args[2]#lung--scale&Anchors;gill--sct$sctAnchors;

###1.data preparation
case_obj <- readRDS("combined_analysis.rds")
DefaultAssay(case_obj) <- "RNA"
case.list <- SplitObject(case_obj, split.by = "orig.ident")
for(i in 1:length(case.list)){
	case.list[[i]]$stim <- "AE_gill"
}
control_obj <- readRDS("combined_analysis.rds")
DefaultAssay(control_obj) <- "RNA"
control.list <- SplitObject(control_obj, split.by = "orig.ident")
for(i in 1:length(control.list)){
        control.list[[i]]$stim <- "FW_gill"
}

Subset_Cells.list <- c(case.list,control.list)

### 2.data integration
if(method == 'lung'){
	for (i in 1:length(Subset_Cells.list)) {
    		Subset_Cells.list[[i]] <- NormalizeData(Subset_Cells.list[[i]], verbose = FALSE)
    		Subset_Cells.list[[i]] <- FindVariableFeatures(Subset_Cells.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
	}
	Integrated <- FindIntegrationAnchors(object.list = Subset_Cells.list, dims = 1:30)
	Integrated <- IntegrateData(anchorset = Integrated, dims = 1:30)
	DefaultAssay(object = Integrated) <- "integrated"
	combined=ScaleData(Integrated,verbose = FALSE)
	combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
}

if(method == 'gill'){
	for (i in 1:length(Subset_Cells.list)) {
        	Subset_Cells.list[[i]] <- SCTransform(Subset_Cells.list[[i]], verbose = FALSE)
        	Subset_Cells.list[[i]] <- NormalizeData(Subset_Cells.list[[i]], normalization.method = "LogNormalize",scale.factor = 10000,assay ="RNA")
        	Subset_Cells.list[[i]] <- ScaleData(Subset_Cells.list[[i]], verbose = FALSE,assay ="RNA")
	}
	obj.features <- SelectIntegrationFeatures(object.list = Subset_Cells.list, nfeatures = 3000)
	obj.list <- PrepSCTIntegration(object.list = Subset_Cells.list, anchor.features = obj.features, verbose = FALSE)

	obj.anchors <- FindIntegrationAnchors(object.list = obj.list, normalization.method = "SCT", anchor.features = obj.features, verbose = FALSE)
	combined <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", verbose = FALSE)
	combined <- RunPCA(object = combined, npcs = 30, verbose = FALSE)
}


###3. normalized

combined <- RunUMAP(object = combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(object = combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
saveRDS(combined,paste0(args[3],".recluster.rds"))
pdf(paste0(args[3],".cluster.pdf"))
DimPlot(object = combined, reduction = "umap",label = TRUE)
dev.off()
par(mfrow=c(3,6))
pdf(paste0(args[3],"dimplot_sample.pdf"),width=24,height=8)
DimPlot(object = combined, reduction = "umap", split.by = "orig.ident")
dev.off()
pdf(paste0(args[3],"dimplot_group.pdf"),width=20,height=8)
DimPlot(object = combined, reduction = "umap", split.by = "stim")
dev.off()
pdf(paste0(args[3],"align_Clusters.pdf"),width=20,height=8)
p1 <- DimPlot(object = combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(object = combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
