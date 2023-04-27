args <- commandArgs(T)
if(length(args) !=2 ){
	cat("\nUsage:\n\t1.cd work_dir\n\t2.run: Rcript 1.FC_mtx.r seurat_obj_rds sample_name\n")
	cat("Note:\n\t1.Packages Require: Seurat\n\t2.Aim: code for calculate FC matrix for specific species\n\n")
	q(save = "no", status = 0, runLast = TRUE)
}

library(Seurat)

Obj_rds <- readRDS(args[1])

###get the FC mt of Object
Objmtx<-as.matrix(Obj_rds@assays$RNA@counts)
ct_geomean=t(apply(Objmtx, 1,  function(x) tapply(x, Obj_rds$cellType,function(y) exp(mean(log(1+y)))-1)))
ct_meansize=tapply(colSums(Objmtx), Obj_rds$cellType, mean)
ideal_cell_size=pmin(1000,median(ct_meansize))
g_fp=t(ideal_cell_size * t(ct_geomean)/as.vector(ct_meansize))
fp_reg=0.05
g_fp_Obj=(fp_reg + g_fp)/apply(fp_reg + g_fp, 1, median)
colnames(g_fp_Obj)<-paste(args[2],colnames(g_fp_Obj),sep="_")
write.table(g_fp_Obj,paste0(args[2],"_broad_FC"),sep="\t",quote=F)

