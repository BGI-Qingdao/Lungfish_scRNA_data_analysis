print(paste("Start time:",format(Sys.time(), "%Y%m%d %X"),sep = " "))

library(Seurat)
library(SoupX)
options(future.globals.maxSize = 100000 * 1024^3)


args<-commandArgs(T)
if(length(args)!=3){
print ("Usage:")
print ("cd work_dir ;Rcript anno.R count_matrix sample max_nFeature_RNA doublet_rate")
print ("Note:")
print ("Aim to run single-sample umap")
q(save = "no", status = 0, runLast = TRUE)
}
setwd(args[1])
#soupx
# toc is flt; tod is raw
toc <- Read10X(args[2],gene.column=1)#RNA_V2
tod <- Read10X(args[3],gene.column=1)## PISA count
tod <- tod[rownames(toc),]

all <- toc
all <- CreateSeuratObject(all)

all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(all)
all <- ScaleData(all, features = all.genes)

all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
all <- FindNeighbors(all, dims = 1:30)
all <- FindClusters(all, resolution = 0.5)
all <- RunUMAP(all, dims = 1:30)

matx <- all@meta.data
sc = SoupChannel(tod, toc)
sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc)
out = adjustCounts(sc,roundToInt=TRUE)
saveRDS(sc,"sc.rds")
DropletUtils:::write10xCounts("soupX_matrix", out,version="3")




