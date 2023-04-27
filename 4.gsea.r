args <- commandArgs(T)
if(length(args) !=5 ){
	cat("\nUsage:\n\tRcript 3.gsea.r seurat_obj_rds MSigDB_ref_path MSigDB_type homologous_gene Out_str\n")
        cat("Note:\n\t1.Env: seurat4\n\t2.Aim: code for gsea enrichment in each cell types\n\n")
        q(save = "no", status = 0, runLast = TRUE)
}

#atacpipline

library(Seurat)
library(clusterProfiler)
library(dplyr)
library(plyr)

object <- readRDS(args[1])
DefaultAssay(object) <- "RNA"
Idents(object) <- object$cell_type1
object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
topn <- object.markers %>% group_by(cluster)
topn <- select(topn, gene, everything())
write.table(topn[topn$p_val_adj <= 0.05,],file=paste0(args[5],"GeneDiffExpFilter.xls"),sep="\t",col.names = TRUE,row.names = F,quote=F)


homo <- read.table(args[4])
names(homo) <- c("gene","symbol")
markers <- object.markers %>% inner_join(homo, by = "gene")

msigdb.GSEA <- list()
msigdb<-read.gmt(args[2])
for(i in 1:length(unique(markers$cluster))){
	k <- markers[which(markers$cluster==unique(markers$cluster)[i]),]
	k1=k[order(k$avg_log2FC,decreasing = T),"avg_log2FC"]
	names(k1) <- k$symbol
	msigdb.GSEA[[i]]=GSEA(k1, TERM2GENE=msigdb, pvalueCutoff = 2,minGSSize=1,maxGSSize=1000)
}
names(msigdb.GSEA) <- as.character(unique(markers$cluster))
df <- ldply(msigdb.GSEA, data.frame)
write.table(df, file = paste0(args[5],"gsea_",args[3],"_result.txt"), sep = "\t", quote = F, row.names = F)
