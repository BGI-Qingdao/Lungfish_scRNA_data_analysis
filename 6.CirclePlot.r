args <- commandArgs(T)
if(length(args) != 10){
        cat("\nUsage:\n\tRcript 2.CirclePlot.r species1_fullname species1_shortname species2_fullname species2_shortname\n\t\t\t      species1_FC_mtx species2_FC_mtx orthologous_pairs ctannPath outdir kldValue\n")
        cat("Note:\n\t1.Conda env: atacpipline\n\t2.Aim: code for circle plot between two species with one to one homologous gene\n\n")
        q(save = "no", status = 0, runLast = TRUE)

}

library(data.table)
library(stringr)
library(ape)
library(tidytree)
library(treeio)
library(dendextend)
library(phangorn)
library(phytools)
library(ComplexHeatmap)
library(eulerr)
library(ggplot2)

scripts <- list.files("~/github/Stylophora_single_cell_atlas-master/metacell_downstream_functions",full.names=TRUE)
for (file in scripts) source(file)

#spfullnames=c(as.character(args[2])=as.character(args[1]),as.character(args[4])=as.character(args[3]))
spfullnames=c(args[1],args[3])
names(spfullnames) <- c(args[2],args[4])
sp1=args[2]
sp2=args[4]
sp=paste(c(sp1,sp2),collapse="_")
sp1_fp_fn=args[5]
sp2_fp_fn=args[6]
OG_pairs_fn=args[7]

csps=csps_create_crossspecies_object(
  sp1_fp_fn = sp1_fp_fn,
  sp2_fp_fn = sp2_fp_fn, 
  OG_pairs_fn = OG_pairs_fn, 
  sp_names = paste0(c(sp1,sp2),"_"), make_sp_colnames = FALSE, quant_norm = TRUE, one2one = TRUE
)
mc_fp=csps$merged
kld <- calcKLD(mc_fp)
kld <- kld$mat
kld<-kld[sort(rownames(kld)),sort(colnames(kld))]

###make color table
#kktable2<-colnames(kld)
#kktable2<-gsub(paste0(args[2],"_"),"",kktable2)
#kktable2<-gsub(paste0(args[4],"_"),"",kktable2)
#kktable<-cbind(colnames(kld),kktable2)
#ctann<-as.data.frame(kktable)
##color1<-rep(c("red","blue"),times=15)
##ctann<-cbind(ctann,color1[1:42])
#sp1length<-length(grep(sp1,colnames(kld)))
#sp2length<-length(colnames(kld))-sp1length
#color1<-c(rainbow(sp1length),rainbow(sp2length))
#ctann<-cbind(ctann,color1)
#colnames(ctann)<-c("metacell","cell_type","color")

ctann <- read.table(args[8], header = T, sep = "\t",comment.char = "")

#Setup circos plotting.  
outdir= args[9]
plotdir <- file.path(outdir,"7.circos")
dir.create(plotdir,showWarnings=FALSE)
sectors.order=ctann$metacell
sector.labels=ctann$cell_type
sectors.order<-as.vector(sectors.order)
sector.labels<-as.vector(sector.labels)

names(sector.labels)=sectors.order
sectors.colors=ctann$color
sectors.colors<-as.vector(sectors.colors)
names(sectors.colors)=sectors.order
sectors.groups=str_extract(sectors.order,paste(sp1,sp2,sep="|"))
names(sectors.groups)=sectors.order

#Plot circos diagram.
circos_file <- file.path(outdir,sprintf("%s_ct_kld.RDS",sp))
suff <- str_remove_all(basename(circos_file),sprintf("%s_|\\.RDS",sp))

threshold=NULL
threshold.quantile=as.numeric(args[10])
annotation.names=TRUE
if(annotation.names==TRUE) {
  sector.groups.padding=c(0.1,0,2.3,0)
  sector.groups.col="gray98"
  sector.groups.border="gray88"
  sector.groups.text.cex=0.68
  sector.groups.text.vjust=-6
  big.gap=5
} else {
  sector.groups.padding=c(0.1,0,0.1,0)
  sector.groups.col="white"
  sector.groups.border="white"
  sector.groups.text.cex=1
  sector.groups.text.vjust=0.5
  big.gap=10
}

name.suffix=sprintf("%s_%s%s",suff,threshold.quantile,ifelse(annotation.names,"","_without_names"))
title.main=suff
title.sub=sprintf("(threshold: %sq)",threshold.quantile)

#
cat("plot begin\n")
chord=csps_plot_chord_diagram(
  mat=kld, sp=c(sp1,sp2), revert=TRUE, 
  regularize.values=FALSE, scale.values=TRUE, threshold=threshold,threshold.quantile=threshold.quantile,
  save.plot=TRUE, save.data=TRUE, sectors.colors=sectors.colors,
  outdir=plotdir, name.suffix=name.suffix, width=9.5, height=9.5,
  sectors.order=sectors.order, sectors.labels=sector.labels,   
  sectors.groups=sectors.groups, sectors.groups.labels=spfullnames, 
  start.degree=0, big.gap=big.gap, annotation.track=TRUE,
  annotation.names=annotation.names, sector.groups.padding=sector.groups.padding,
  self.links=FALSE, sector.text.cex=0.38, 
  sector.groups.col=sector.groups.col, sector.groups.border=sector.groups.border,
  sector.groups.text.cex=sector.groups.text.cex, 
  sector.groups.text.vjust=sector.groups.text.vjust
  #title.main=title.main, title.sub=title.sub
)
circos.clear()

saveRDS(mc_fp,file.path(plotdir,sprintf("circos_%s_ct_csps.RDS",sp)))

