args <- commandArgs(T)
if(length(args) != 3){
	cat("\nUsage:\n\tRscript 1.data_pre.r seurat_obj_rds eggnog_orth_list obj_stim\n")
	cat("Note:\n\t1.Packages Require:Seurat\n\t2.Aim:code for prepare matrix to run cellphonedb\n")
	q(save = "no", status = 0, runLast = TRUE)
}

library(Seurat)

sub_feature_seurat <- function(homogene, seurat_obj){
        counts <- GetAssayData(seurat_obj, assay = "RNA")
        counts <- counts[(which(rownames(counts) %in% homogene$V1)),]
        sub_obj <- subset(seurat_obj, features = rownames(counts))
        geneName <- rownames(sub_obj@assays$RNA@counts)
        for(i in 1:length(geneName)){if(geneName[i] %in% homogene$V1) geneName[i]=homogene[which(homogene$V1==geneName[i]),2]}
        rownames(sub_obj@assays$RNA@counts) <- geneName
	rownames(sub_obj@assays$RNA@data) <- geneName
        return(sub_obj)
}

object <- readRDS(args[1])
homoGene <- read.table(args[2])
DefaultAssay(object) <- "RNA"
Idents(object) <- object$stim
sub_obj <- subset(object, idents = args[3])
sub_obj <- NormalizeData(sub_obj, normalization.method = "LogNormalize", scale.factor = 10000)
sub_obj <- sub_feature_seurat(homoGene, sub_obj)
write.table(as.matrix(sub_obj@assays$RNA@data), 'cellphonedb_count.txt', sep='\t', quote=F)

meta_data <- cbind(rownames(sub_obj@meta.data), sub_obj@meta.data[,'cell_type1', drop=F])
names(meta_data) <- c("Cell", "cell_type")
meta_data <- as.matrix(meta_data)
meta_data[is.na(meta_data)] = "Unknown" #There can't be... In the cell type NA
write.table(meta_data, 'cellphonedb_meta.txt', sep='\t', quote=F, row.names=F)
