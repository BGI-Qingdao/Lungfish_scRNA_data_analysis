args <- commandArgs(T)
if(length(args)!=5){
	cat("\nUsage:\n\tRscript net_plot.r count_network_file count_threshold lable_cex lable_dist prefix\n")
	cat("Note:\n\t1.conda_env:seurat4\n\t2.Aim:code for cellphonedb visualize\n")
	q(save = "no", status = 0, runLast = TRUE)
}

library(psych)
library(qgraph)
library(igraph)
library(purrr)
library(ggplot2)

mynet <- read.delim(args[1], check.names = FALSE)
mynet <- mynet[which(mynet$count > as.numeric(args[2])),]
vlcex = as.numeric(args[3]) 
vldist = as.numeric(args[4]) 
net<- graph_from_data_frame(mynet)
V(net)$size <- Matrix::rowSums(net[,])*2
karate_groups <- cluster_optimal(net)
coords <- layout_in_circle(net, order = order(membership(karate_groups)))

pal <- c("#882E72", "#B178A6", "#D4A6C8", "#FABFD2", "#c8262e", "#ff6347",
         "#f5a9a4", "#4EB265", "#90C987", "#CAE0AB", "#405AB8", "#7AA1DD",
         "#536872")
names(pal) <- c("Alveolar cell",
                 "Epithelial cell",
                 "Goblet cell",
                 "Ciliated cell",
                 "Erythroid cell",
                 "Vasc_Endo",
                 "Vasc_SMC",
                 "Stromal cell",
                 "Stroaml cell_CLDN",
                 "SMC",
                 "Lymphoid cell",
                 "Macrophage",
                 "Immune cell_F13A1"
                 )
pdf(paste0(args[5],"_net.pdf"), width=24, height=16)
#plot
#----------------------------
V(net)$color <- pal[names(V(net))]
E(net)$width  <- E(net)$count*3
net2 <- net

edge.start <- ends(net, es=E(net), names=F)[,1]
edge.col <- V(net)$color[edge.start]

p <- plot(net, edge.arrow.size=.1,
     	edge.curved=0.1,
     	vertex.size = 5*Matrix::rowSums(net[,]),
     	vertex.frame.color="#555555",
     	vertex.label.color="black",
     	vertex.label.cex = vlcex,
     	vertex.label.dist = vldist,
     	layout = coords,edge.color=edge.col)

print(p)

dev.off()

