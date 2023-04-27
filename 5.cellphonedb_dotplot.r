args <- commandArgs(T)
if(length(args) != 6){
	cat("\nUsage:\n\tRscript bubble_plot.r pval_file mean_file selected_celltype selected_LRpairs filter_mean prefix \n")
	cat("Notion:\n\t1.conda_env:seurat4\n\t2.aim:code for plot bubble diagram with cellphonedb result")
	q(save = "no", status = 0, runLast = TRUE)
}

library(psych)
library(qgraph)
library(igraph)
library(ggplot2)

mypvals <- read.delim(args[1], check.names = FALSE)
mymeans <- read.delim(args[2], check.names = FALSE)
cell_type = args[3]
LR_pairs = args[4]
if(cell_type == "all"){
	kp = as.logical(rep("TRUE",ncol(mypvals)))
}else{
	kp = grepl(pattern = cell_type, colnames(mypvals))
}
#table(kp)
pos = (1:ncol(mypvals))[kp] 
choose_pvalues <- mypvals[,c(c(1,5,6,8,9),pos)]
choose_means <- mymeans[,c(c(1,5,6,8,9),pos)]
logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<0.05, 1, sum) 

choose_pvalues <- choose_pvalues[logi>=1,]

logi1 <- choose_pvalues$gene_a != ""
logi2 <- choose_pvalues$gene_b != ""
logi <- logi1 & logi2
choose_pvalues <- choose_pvalues[logi,]


choose_means <- choose_means[choose_means$id_cp_interaction %in% choose_pvalues$id_cp_interaction,]

dim(choose_means)
dim(choose_pvalues)


library(tidyverse)
meansdf <- choose_means %>% reshape2::melt()
meansdf <- data.frame(interacting_pair = paste0(meansdf$gene_a,"_",meansdf$gene_b),
                      CC = meansdf$variable,
                      means = meansdf$value)
pvalsdf <- choose_pvalues %>% reshape2::melt()
pvalsdf <- data.frame(interacting_pair = paste0(pvalsdf$gene_a,"_",pvalsdf$gene_b),
                      CC = pvalsdf$variable,
                      pvals = pvalsdf$value)

pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
pldf <- merge(pvalsdf,meansdf,by = "joinlab")

if(LR_pairs == "all"){
}else{
	select_pairs <- grep(LR_pairs, pldf$interacting_pair.x, value = T)
	pldf%>% dplyr::filter(interacting_pair.x %in% select_pairs) -> pldf
}


# dotplot
summary((dplyr::filter(pldf,means >0))$means)
#head(pldf)
pcc = pldf%>% dplyr::filter(means > as.numeric(args[5])) %>%  dplyr::filter(pvals < 0.05) %>%
  ggplot(aes(CC.x,interacting_pair.x) )+ 
  geom_point(aes(color=means,size=-log10(pvals+0.0001)) ) +
  scale_size_continuous(range = c(1,3))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue",midpoint = 2  )+ 
  theme_bw()+ 
  # scale_color_manual(values = rainbow(100))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))

pdf(paste0(args[6], "_mean_", args[5], "_bubble.pdf"), width=20, height=12)
print(pcc)
dev.off()

write.table(pldf, file = paste0(args[6], "_bubble.txt"), sep = "\t", quote = F, row.names = F)

